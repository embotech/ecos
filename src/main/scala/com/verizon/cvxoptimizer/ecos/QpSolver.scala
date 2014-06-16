package com.verizon.cvxoptimizer.ecos

import breeze.linalg.DenseMatrix
import breeze.linalg.CSCMatrix
import breeze.linalg.max
import breeze.optimize.LBFGS
import breeze.linalg.DenseVector
import breeze.optimize.DiffFunction
import breeze.linalg.norm
import org.jblas.Solve
import org.jblas.DoubleMatrix
import org.jblas.Decompose

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

/*
 * Send the dimension of Hessian n, Aeq (structure of equality constraints) and A (structure of inequality constraints)
 */
class QpSolver(n: Int, diagonal: Boolean = false,
  Equalities: Option[CSCMatrix[Double]] = None, Inequalities: Option[CSCMatrix[Double]] = None,
  lbFlag: Boolean = false, ubFlag: Boolean = false) {
  NativeECOS.loadLibraryAndCheckErrors()
  
  //Append stacked variables a, b and t to x to generate SOCP variables
  val c = Array.fill[Double](n + 3)(0.0)
  c.update(n, 0)
  c.update(n + 1, 0)
  c.update(n + 2, 1)

  var bounds = 0;
  if (lbFlag) bounds += 1;
  if (ubFlag) bounds += 1;
  
  //Aeqx = Beq is represented by linearEquality x n
  val linearEquality = if (Equalities == None) { 0 } else Equalities.get.rows

  val aeqBuilder = new CSCMatrix.Builder[Double](linearEquality + 2, n + 3)
  aeqBuilder.add(0, n, 1)
  aeqBuilder.add(0, n + 2, -0.5)
  aeqBuilder.add(1, n + 1, 1)
  aeqBuilder.add(1, n + 2, -0.5)

  if (linearEquality > 0) {
    for (i <- Equalities.get.iterator) {
      val ((row, col), value) = i
      aeqBuilder.add(2 + row, col, value)
    }
  }

  val beqBuilder = Array.fill[Double](linearEquality + 2)(0.0)
  beqBuilder.update(linearEquality, -0.5)
  beqBuilder.update(linearEquality + 1, 0.5)

  //Ax <= B is represented by linearInequality x n
  val linearInequality = if (Inequalities == None) { 0 } else Inequalities.get.rows

  //Ax <= B is represented by linearInequality x n
  //upper bound is represented by n x n sparse identity
  //lower bound is represented by n x n sparse identity

  val gBuilder = new CSCMatrix.Builder[Double](linearInequality + bounds * n + n + 2, n + 3)

  if (ubFlag) {
    for (j <- 0 to n - 1) gBuilder.add(linearInequality + j, j, 1)
  }

  if (lbFlag) {
    if (ubFlag) for (j <- 0 to n - 1) gBuilder.add(linearInequality + n + j, j, -1)
    else for (j <- 0 to n - 1) gBuilder.add(linearInequality + j, j, -1)
  }

  //Create workspace for G
  gBuilder.add(linearInequality + bounds * n, n + 1, -math.sqrt(2))
  if(diagonal) {
    for(diag <- 0 to n-1) gBuilder.add(linearInequality + bounds*n + 1 + diag, diag, 0)
  } 
  else {
    for (row <- 0 to n - 1)
      for (col <- 0 to n - 1)
      gBuilder.add(linearInequality + bounds * n + 1 + row, col, 0)
  }
  gBuilder.add(linearInequality + bounds * n + n + 1, n, -math.sqrt(2))

  val hBuilder = Array.fill[Double](linearInequality + bounds * n + n + 2)(0.0)

  if (linearInequality > 0) {
    for (i <- Inequalities.get.iterator) {
      val ((row, col), value) = i
      gBuilder.add(row, col, value)
    }
  }

  val G = gBuilder.result
  val Aeq = aeqBuilder.result

  val (linear, cones) = if (linearInequality > 0) {
    (linearInequality, Array[Int](n + 2))
  } else {
    (0, Array[Int](n + 2))
  }

  //Used for recomputing the solution
  var x = Array.fill[Double](c.length)(0.0)

  def updateEquality(beq: Array[Double]) {
    if (linearEquality != beq.length)
      throw new IllegalArgumentException("QpSolver: mismatch on equality constraints")
    for (i <- 0 to beq.length - 1) beqBuilder(2 + i) = beq(i)
  }

  def updateInequality(A: CSCMatrix[Double], b: Array[Double]) {
    if (linearInequality != b.length) {
      throw new IllegalArgumentException("QpSolver: mismatch on inequality constraints")
    }
    for (i <- 0 to b.length - 1) hBuilder(i) = b(i)
  }
  
  def updateDiagonal(diag: Array[Double]) {
    val w = diag.map(d=>math.sqrt(d))
    for(d <- 0 to w.length-1) G.update(linearInequality + bounds * n + 1 + d, d, -w(d))
  }
  
  def updateHessian(H: DoubleMatrix) {
    val W = Decompose.cholesky(H)
    if (n != W.rows) {
      throw new IllegalArgumentException("QpSolver Hessian rows should be same as workspace rows")
    }
    for (row <- 0 to n - 1)
      for (col <- 0 to n - 1) {
        val w = W.get(row, col)
        G.update(linearInequality + bounds * n + 1 + row, col, -w)
      }
  }

  //TODO : Should we add any asserts here or it should be clean
  def updateCholesky(row: Int, col: Int, value: Double) {
    require(row < n)
    require(col < n)
    G.update(linearInequality + bounds * n + 1 + row, col, -value)
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
  
  def run(H: DoubleMatrix, f: Array[Double]): (Int, Array[Double]) = {
    if(diagonal) {
      throw new IllegalArgumentException("Qpsolver: digonal flag must be false for dense solve")
    }
    val hessianStart = System.currentTimeMillis()
    updateHessian(H)
    val hessianTime = System.currentTimeMillis() - hessianStart
    
    val linearStart = System.currentTimeMillis()
    updateLinearObjective(f)
    val linearTime = System.currentTimeMillis() - linearStart
    
    val nativeStart = System.currentTimeMillis()
    val status = NativeECOS.solveSocp(c, G.rows, G.cols, G.data, G.colPtrs, G.rowIndices, hBuilder,
      Aeq.rows, Aeq.cols, Aeq.data, Aeq.colPtrs, Aeq.rowIndices, beqBuilder,
      linear, cones, x);
    val nativeTime = System.currentTimeMillis() - nativeStart
    println("hessian " + hessianTime + " linear " + linearTime + " native " + nativeTime)
    
    (status, x.slice(0, n))
  }
  
  def run(Hdiag: Array[Double], f: Array[Double]) : (Int, Array[Double]) = {
    val diagStart = System.currentTimeMillis()
    if (!diagonal) {
      throw new IllegalArgumentException("QpSolver: diagonal flag must be true for sparse solve")
    }
    updateDiagonal(Hdiag)
    val diagTime = System.currentTimeMillis() - diagStart

    val linearStart = System.currentTimeMillis()
    updateLinearObjective(f)
    val linearTime = System.currentTimeMillis() - linearStart

    val nativeStart = System.currentTimeMillis()
    val status = NativeECOS.solveSocp(c, G.rows, G.cols, G.data, G.colPtrs, G.rowIndices, hBuilder,
      Aeq.rows, Aeq.cols, Aeq.data, Aeq.colPtrs, Aeq.rowIndices, beqBuilder,
      linear, cones, x);
    val nativeTime = System.currentTimeMillis() - nativeStart
    println("diagonal " + diagTime + " linear " + linearTime + " native " + nativeTime)
    
    (status, x.slice(0, n - 1))
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
    println("Test jecos QpSolver, breeze LBFGS and jblas solvePositive with problemSize " + problemSize)

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
    
    val Hdiag = Array.fill[Double](problemSize)(2.0)
    val H = DoubleMatrix.eye(problemSize).mul(2.0)
    val f = Array.fill[Double](problemSize)(-6)

    val qpSolverDiag = new QpSolver(problemSize, true)
    val qpDiagStart = System.currentTimeMillis()
    val (statusDiag, qpResultDiag) = qpSolverDiag.run(Hdiag, f)
    val qpDiagTime = System.currentTimeMillis() - qpDiagStart
    
    val qpSolver = new QpSolver(problemSize)
    val qpStart = System.currentTimeMillis()
    val (status, qpResult) = qpSolver.run(H, f)
    val qpTime = System.currentTimeMillis() - qpStart
    
    println("Qp output")
    for (i <- 0 to qpResult.length - 1) println(qpResult(i))

    val jblasH = DoubleMatrix.eye(problemSize).mul(2.0)
    val jblasf = DoubleMatrix.zeros(problemSize, 1).add(6.0)

    val dposvStart = System.currentTimeMillis()
    val dposvResult = Solve.solvePositive(jblasH, jblasf).data
    println("Dposv output")
    for (i <- 0 to dposvResult.length - 1) println(dposvResult(i))
    
    val dposvTime = System.currentTimeMillis() - dposvStart

    println("Runtime bfgs " + bfgsTime + " qpDiag " + qpDiagTime + " qp " + qpTime + " dposv " + dposvTime)
    
    val n = 5
    val ata = new DoubleMatrix(Array(
      Array( 4.377, -3.531, -1.306, -0.139,  3.418),
      Array(-3.531,  4.344,  0.934,  0.305, -2.140),
      Array(-1.306,  0.934,  2.644, -0.203, -0.170),
      Array(-0.139,  0.305, -0.203,  5.883,  1.428),
      Array( 3.418, -2.140, -0.170,  1.428,  4.684)))
    val atb = new DoubleMatrix(Array(-1.632, 2.115, 1.094, -1.025, -0.636))

    val goodx = Array(0.13025, 0.54506, 0.2874, 0.0, 0.028628)
    
    val qpSolverBounds = new QpSolver(n, false, None, None, true, false)
    
    val (statusBounds, qpResultBounds) = qpSolverBounds.run(ata, atb.data)
    
    for(i <- 0 until n) {
      println(qpResultBounds(i) + " " + goodx(i))
      //assert(Math.abs(x(i) - goodx(i)) < 1e-3)
      //assert(x(i) >= 0)
    }
  }
}