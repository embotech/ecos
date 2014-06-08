package com.verizon.jecos

import breeze.linalg.DenseMatrix
import breeze.linalg.CSCMatrix
import breeze.linalg.max
import breeze.linalg.cholesky
import breeze.optimize.LBFGS
import breeze.linalg.DenseVector
import breeze.optimize.DiffFunction
import breeze.linalg.norm
import org.jblas.Solve
import org.jblas.DoubleMatrix

/* QpSolver solves quadratic programming problem as follows:
 *   
 * min 0.5*x'*H*x + f'*x 
 * subject to: A*x <= b 
 * while additionally satisfying the equality constraint Aeq*x = beq and defines a set of lower and upper
 * bounds on the variables x so that lb <= x <= ub
 * 
 * If no bounds exist use empty matrices for LB and UB Set LB(i) = -Inf if
 * X(i) is unbounded below; Set UB(i) = Inf if X(i) is unbounded above.
 * 
 * QpSolver rewrites the QP as a second order cone program (SOCP) that can
 * be solved by ECOS. 
 * 
 * The rewrite introduces 3 variables t, a, b 
 * 	 min f'x + t
 * x,t,a,b 
 * 	 s.t. A*x <= b 
 *   	[ lb <= x <= ub ] 
 *    	[ Aeq*x == beq ] 
 *     	c = (t-1)/2   }
 * 	   	d = (t+1)/2   }  these constraints say 0.5*x'*H*x <= t
 *      ||W*x||       }  where W'*W = 0.5*H
 *      || c ||2 <= d }
 *      ||W*x || } || c ||2 <=
 *
 * For unconstrained quadratic program (least squares) an interesting comparison will be with the Ceres solver from Google
 * http://ceres-solver.org

 * 
 * min 0.5*x'*H*x + f'*x 
 * s.t 
 * 		Ax <= B (to be used for epigraph formulation of L1 constraint)
 * 		AEQx = Beq (To be used in Orthogonality / Orthonormality constraint)
 * 		[lb <= x <= ub] optional
 * 
 * See jecos/matlab/ecosqp.m for further details
 * 
 */
class QpSolver(val H: DenseMatrix[Double], val f: Array[Double],
  val linearInequality: Int, val linearEquality: Int,
  val lbFlag: Boolean, val ubFlag: Boolean) {

  NativeECOS.loadLibraryAndCheckErrors()

  //TO DO : The cholesky factorization default returns the L matrix, breeze has to be updated to U matrix 
  val W = cholesky(H)
  
  val n = H.rows
  
  if (n != f.length) {
    throw new IllegalArgumentException("QpSolver gram matrix rows should be equal to linear objective dimension")
  }
  
  //Append stacked variables a, b and t to x to generate SOCP variables
  val c = f ++ Array[Double](0, 0, 1)
  
  var bounds = 0;
  if (lbFlag) bounds += 1;
  if (ubFlag) bounds += 1;

  println("Lower bound " + lbFlag + " Upper bound " + ubFlag)

  //Aeqx = Beq is represented by linearEquality x n
  val aeqBuilder = new CSCMatrix.Builder[Double](linearEquality + 2, n + 3)
  aeqBuilder.add(0, n, 1)
  aeqBuilder.add(0, n + 2, -0.5)
  aeqBuilder.add(1, n + 1, 1)
  aeqBuilder.add(1, n + 2, -0.5)
  
  val beqConst = Array[Double](-0.5, 0.5)
  val beqBuilder = if(linearEquality > 0 )  beqConst ++ Array[Double](linearEquality) 
  				   else beqConst
  
  //Ax <= B is represented by linearInequality x n
  //upper bound is represented by n x n sparse identity
  //lower bound is represented by n x n sparse identity
  /*
  val aBuilder = new CSCMatrix.Builder[Double](linearInequality + bounds * n, n)

  if (ubFlag) {
    for (j <- 0 to n - 1) aBuilder.add(linearInequality + j, j, 1)
  }

  if (lbFlag) {
    if (ubFlag) for (j <- 0 to n - 1) aBuilder.add(linearInequality + n + j, j, -1)
    else for (j <- 0 to n - 1) aBuilder.add(linearInequality + j, j, -1)
  }

  val bBuilder = new Array[Double](linearInequality + bounds * n)

  if (ubFlag) {
    for (j <- 0 to n - 1) bBuilder(linearInequality + j)
  }
  * 
  */
  val gBuilder = new CSCMatrix.Builder[Double](linearInequality + bounds * n + n + 2, n + 3)
  
  if (ubFlag) {
    for (j <- 0 to n - 1) gBuilder.add(linearInequality + j, j, 1)
  }

  if (lbFlag) {
    if (ubFlag) for (j <- 0 to n - 1) gBuilder.add(linearInequality + n + j, j, -1)
    else for (j <- 0 to n - 1) gBuilder.add(linearInequality + j, j, -1)
  }

  gBuilder.add(linearInequality + bounds*n, n + 1, -math.sqrt(2))
  for(i <- W.iterator) {
    val ((row, col), value) = i
    if ( value != 0) gBuilder.add(linearInequality + bounds*n + 1 + row, col, -value) 
  }
  gBuilder.add(linearInequality + bounds*n + n + 1, n, -math.sqrt(2))
  
  val hBuilder = Array.fill[Double](linearInequality + bounds * n + n + 2)(0)
  
  def updateEquality(Aeq: CSCMatrix[Double], beq: Array[Double]) {
    if (linearEquality != beq.length)
      throw new IllegalArgumentException("QpSolver: mismatch on equality constraints")
    val iter = Aeq.activeIterator
    iter.map { case ((row, col), value) => aeqBuilder.add(2 + row, col, value) }
    for (i <- 0 to beq.length - 1) beqBuilder(2 + i) = beq(i)
  }
  def updateInequality(A: CSCMatrix[Double], b: Array[Double]) {
    if (linearInequality != b.length) {
      throw new IllegalArgumentException("QpSolver: mismatch on inequality constraints")
    }
    val iter = A.activeIterator
    iter.map { case ((row, col), value) => gBuilder.add(row, col, value) }
    for (i <- 0 to b.length - 1) hBuilder(i) = b(i)
  }
  def updateLb(lb: Array[Double]) {
    require(lbFlag)
    if (ubFlag) for (j <- 0 to n - 1) hBuilder(linearInequality + n + j) = -lb(j)
    else for (j <- 0 to n - 1) hBuilder(linearInequality + j) = lb(j)
  }
  def updateUb(ub: Array[Double]) {
    require(ubFlag)
    for (j <- 0 to n - 1) hBuilder(linearInequality + j) = ub(j)
  }
  def updateLinearObjective(f: Array[Double]) {
    if (f.length != c.length - 3) {
      throw new IllegalArgumentException("QpSolver: mismatch on dimension on linear objective update")
    }
    for (i <- 0 to f.length - 1) c.update(i, f(i))
  }
  def run: (Int, Array[Double]) = {
    val G = gBuilder.result
    val Aeq = aeqBuilder.result
    
    val (linear, cones) = if (linearInequality > 0) {
      (linearInequality, Array[Int](n + 2))
    } else {
      (0, Array[Int](n + 2))
    }
    
    val x = Array.fill[Double](c.length)(0.0)
    
    val status = NativeECOS.solveSocp(c, G.rows, G.cols, G.data, G.colPtrs, G.rowIndices, hBuilder,
      Aeq.rows, Aeq.cols, Aeq.data, Aeq.colPtrs, Aeq.rowIndices, beqBuilder,
      linear, cones, x);
    
    return (status, x)
  }
}

object QpSolver {
  def main(args: Array[String]) {    
    if (args.length < 1) {
      println("Usage: QpSolver 100")
      println("Test QpSolver with a simple quadratic function of dimension 100")
      sys.exit(1)
    }
    
    val problemSize = args(0).toInt
    println("Test jecos QpSolver, breeze LBFGS and jblas solvePositive with problemSize" + problemSize)

    val lbfgs = new LBFGS[DenseVector[Double]](problemSize, 4)
    def optimizeThis(init: DenseVector[Double]) = {
      val f = new DiffFunction[DenseVector[Double]] {
        def calculate(x: DenseVector[Double]) = {
          (norm((x - 3.0) :^ 2.0, 1), (x * 2.0) - 6.0)
        }
      }
      val result = lbfgs.minimize(f, init)
      norm(result - 3.0, 2) < 1E-10
    }
    
    val init = DenseVector.zeros[Double](problemSize)
    init(0 until init.length by 2) := -1.2
    init(1 until init.length by 2) := 1.0
    
    val startBfgs = System.currentTimeMillis()
    assert(optimizeThis(init))
    val bfgsTime = System.currentTimeMillis() - startBfgs
    
    val H = DenseMatrix.eye[Double](problemSize)*2.0    
    val f = Array.fill[Double](problemSize)(-6)
    
    val qpSolver = new QpSolver(H, f, 0, 0, false, false)
    val qpStart = System.currentTimeMillis()
    val (status, qpResult) = qpSolver.run
    val qpTime = System.currentTimeMillis() - qpStart
    println("Qp output")
    for(i <- 0 to qpResult.length-1) println(qpResult(i))
    
    val jblasH = DoubleMatrix.eye(problemSize).mul(2.0)    
    val jblasf = DoubleMatrix.zeros(problemSize, 1).add(6.0)
    
    val dposvStart = System.currentTimeMillis()
    val dposvResult = Solve.solvePositive(jblasH, jblasf).data
    println("Dposv output")
    for(i <- 0 to dposvResult.length - 1) println(dposvResult(i))
    
    val dposvTime = System.currentTimeMillis() - dposvStart
    
    println("Runtime bfgs " + bfgsTime + " qp " + qpTime + " dposv " + dposvTime)
  }
}