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
import scala.math.sqrt
import scala.math.abs
import breeze.optimize.OWLQN

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
 * The rewrite introduces a variable t such that
 * 	 min t
 *   x,t
 * 	 s.t. A*x <= b 
 *   	[ lb <= x <= ub ] 
 *    	[ Aeq*x == beq ] 
 *     ||         Wx         ||
 *     || (t-f'*x+1)/sqrt(2) ||_2 <= (t - f'*x - 1)/sqrt(2)
 *     
 * For unconstrained quadratic program (least squares) an interesting comparison will be with the Ceres solver from Google
 * http://ceres-solver.org. Look at how they are handling dense Hessian ? 
 * 
 * See cvxoptimizer/matlab/ecosqp.m for further details
 * 
 */

/*
 * Send the dimension of Hessian n, Aeq (structure of equality constraints) and A (structure of inequality constraints)
 * 
 * nHessian : variables which are part of Hessian
 * 
 * nLinear : variables which are part of linear objective but not hessian
 */
class QpSolver(nHessian: Int, nLinear: Int = 0, diagonal: Boolean = false,
  Equalities: Option[CSCMatrix[Double]] = None, Inequalities: Option[CSCMatrix[Double]] = None,
  lbFlag: Boolean = false, ubFlag: Boolean = false) {
  NativeECOS.loadLibraryAndCheckErrors()
  
  val n = nHessian + nLinear
  
  //The new variable is stacked as [x, t]
  val c = Array.fill[Double](n + 1)(0.0)
  c.update(n, 1.0)
  
  var bounds = 0;
  if (lbFlag) bounds += 1;
  if (ubFlag) bounds += 1;
  
  //Aeqx = Beq is represented by linearEquality x n
  val linearEquality = if (Equalities == None) { 0 } else Equalities.get.rows

  //pad Aeq with a zero column for t
  val aeqBuilder = new CSCMatrix.Builder[Double](linearEquality, n + 1)

  if (linearEquality > 0) {
    for (i <- Equalities.get.iterator) {
      val ((row, col), value) = i
      aeqBuilder.add(row, col, value)
    }
  }
  
  val beqBuilder = Array.fill[Double](linearEquality)(0.0)
  
  //Ax <= B is represented by linearInequality x n
  val linearInequality = if (Inequalities == None) 0 else Inequalities.get.rows
  
  //Ax <= B is represented by linearInequality x n
  //upper bound is represented by n x n sparse matrix
  //lower bound is represented by n x n sparse matrix

  val gBuilder = new CSCMatrix.Builder[Double](linearInequality + bounds*n + n + 2, n + 1)
  if (ubFlag) for (j <- 0 until n) gBuilder.add(linearInequality + j, j, 1)

  if (lbFlag) {
    if (ubFlag) for (j <- 0 until n) gBuilder.add(linearInequality + n + j, j, -1)
    else for (j <- 0 until n) gBuilder.add(linearInequality + j, j, -1)
  }
  
  //Create workspace for G
  /*
  	Gquad = [fhalf', -1/sqrt(2);
         	-W, zerocolumn;
         	-fhalf', +1/sqrt(2)];
	hquad = [1/sqrt(2); zerocolumn; 1/sqrt(2)];
  *
  */
  var i = 0
  while (i < n) {
    gBuilder.add(linearInequality + bounds * n, i, 0)
    i = i + 1
  }
  gBuilder.add(linearInequality + bounds * n, n, -1 / sqrt(2))

  if (diagonal) {
    i = 0
    while (i < n) {
      gBuilder.add(linearInequality + bounds * n + 1 + i, i, 0.0)
      i = i + 1
    }
  } else {
    for (row <- 0 to n - 1)
      for (col <- 0 to n - 1)
        gBuilder.add(linearInequality + bounds * n + 1 + row, col, 0.0)
  }

  i = 0
  while (i < n) {
    gBuilder.add(linearInequality + bounds * n + n + 1, i, 0)
    i = i + 1
  }
  gBuilder.add(linearInequality + bounds * n + n + 1, n, 1 / sqrt(2))

  /* Default lower and upper bounds are 0.0, they have to be set explicitly */
  val hBuilder = Array.fill[Double](linearInequality + bounds * n + n + 2)(0.0)
  
  hBuilder.update(linearInequality + bounds * n, 1.0/sqrt(2))
  hBuilder.update(linearInequality + bounds * n + n + 1, 1.0/sqrt(2))

  if (linearInequality > 0) {
    println("Update inequality")
    for (i <- Inequalities.get.iterator) {
      val ((row, col), value) = i
      gBuilder.add(row, col, value)
    }
  }
  
  val G = gBuilder.result
  val Aeq = aeqBuilder.result
  
  val (linear, cones) = if (linearInequality > 0) {
    (linearInequality + bounds*n, Array[Int](n + 2))
  } else {
    (bounds*n, Array[Int](n + 2))
  }
  
  //Used for recomputing the solution
  var x = Array.fill[Double](c.length)(0.0)

  def updateEquality(beq: Array[Double]) {
    if (linearEquality != beq.length)
      throw new IllegalArgumentException("QpSolver: mismatch on equality constraints")
    for (i <- 0 to beq.length - 1) beqBuilder.update(i,beq(i))
  }
  
  def updateInequality(A: CSCMatrix[Double], b: Array[Double]) {
    if (linearInequality != b.length) {
      throw new IllegalArgumentException("QpSolver: mismatch on inequality constraints")
    }
    for (i <- 0 to b.length - 1) hBuilder.update(i,b(i))
  }
  
  def updateDiagonal(diag: Array[Double]) {
    val w = diag.map(d => math.sqrt(d))
    for (d <- 0 to w.length - 1) G.update(linearInequality + bounds * n + 1 + d, d, -w(d))
  }

  def updateHessian(H: DoubleMatrix) {
    val W = Decompose.cholesky(H)
    if (n != W.rows) {
      throw new IllegalArgumentException("QpSolver Hessian rows should be same as workspace rows")
    }
    for (row <- 0 to nHessian - 1)
      for (col <- 0 to nHessian - 1) {
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
    if (!lbFlag) {
      throw new IllegalArgumentException("QpSolver use updateLb if solver created with lower bound feature")
    }
    if (lb.length != n) {
      throw new IllegalArgumentException("QpSolver problem size should be same as lower bound length")
    }
    if (ubFlag) for (j <- 0 to n - 1) hBuilder.update(linearInequality + n + j, -lb(j))
    else for (j <- 0 to n - 1) hBuilder.update(linearInequality + j, -lb(j))
  }

  def updateUb(ub: Array[Double]) {
    if(!ubFlag) {
      throw new IllegalArgumentException("QpSolver use updateUb if solver created with upper bound feature")
    }
    if (ub.length != n) {
      throw new IllegalArgumentException("QpSolver problem size should be same as upper bound length")
    }
    for (j <- 0 to n - 1) hBuilder.update(linearInequality + j, ub(j))
  }

  def updateLinearObjective(f: Array[Double]) {
    if (f.length != c.length - 1) {
      throw new IllegalArgumentException("QpSolver: mismatch on dimension on linear objective update")
    }
    val ftrans = f.map(x => x / sqrt(2))

    //Update G
    for(i <- 0 until f.length) {
      G.update(linearInequality + bounds * n, i, ftrans(i))
      G.update(linearInequality + bounds * n + n + 1, i, -ftrans(i))
    }
  }

  def run(H: DoubleMatrix, f: Array[Double]): (Int, Array[Double]) = {
    if (diagonal) {
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

  def run(Hdiag: Array[Double], f: Array[Double]): (Int, Array[Double]) = {
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
    
    (status, x.slice(0, nHessian))
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
    
    def optimizeWithLBFGS(init: DenseVector[Double]) = {
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
    assert(optimizeWithLBFGS(init))
    val bfgsTime = System.currentTimeMillis() - startBfgs

    val jblasH = DoubleMatrix.eye(problemSize).mul(2.0)
    val jblasf = DoubleMatrix.zeros(problemSize, 1).add(6.0)

    val dposvStart = System.currentTimeMillis()
    val dposvResult = Solve.solvePositive(jblasH, jblasf).data
    val dposvTime = System.currentTimeMillis() - dposvStart

    val Hdiag = Array.fill[Double](problemSize)(2.0)
    val H = DoubleMatrix.eye(problemSize).mul(2.0)
    val f = Array.fill[Double](problemSize)(-6)

    val qpSolverDiag = new QpSolver(problemSize, 0, true)
    val qpDiagStart = System.currentTimeMillis()
    val (statusDiag, qpResultDiag) = qpSolverDiag.run(Hdiag, f)
    val qpDiagTime = System.currentTimeMillis() - qpDiagStart

    val qpSolver = new QpSolver(problemSize)
    val qpStart = System.currentTimeMillis()
    val (status, qpResult) = qpSolver.run(H, f)
    val qpTime = System.currentTimeMillis() - qpStart
    
    //TO DO : Except basic runtime comparisons move everything else to tests
    assert(norm(DenseVector(qpResult) - DenseVector(dposvResult), 2) < 1E-3)
    assert(norm(DenseVector(qpResultDiag) - DenseVector(dposvResult), 2) < 1E-3)

    println("Runtime bfgs " + bfgsTime + " qpDiag " + qpDiagTime + " qp " + qpTime + " dposv " + dposvTime)

    val n = 5
    val ata = new DoubleMatrix(Array(
      Array(4.377, -3.531, -1.306, -0.139, 3.418),
      Array(-3.531, 4.344, 0.934, 0.305, -2.140),
      Array(-1.306, 0.934, 2.644, -0.203, -0.170),
      Array(-0.139, 0.305, -0.203, 5.883, 1.428),
      Array(3.418, -2.140, -0.170, 1.428, 4.684)))
    
    val atb = new DoubleMatrix(Array(-1.632, 2.115, 1.094, -1.025, -0.636))
    
    val goodx = Array(0.13025, 0.54506, 0.2874, 0.0, 0.028628)
    
    /*QpSolver with lower bound true, upper bound false*/
    val qpSolverBounds = new QpSolver(n, 0, false, None, None, true, false)
    
    val (statusBounds, qpResultBounds) = qpSolverBounds.run(ata, atb.mul(-1).data)
    
    for (i <- 0 until n) {
      println(qpResultBounds(i) + " " + goodx(i))
      assert(Math.abs(qpResultBounds(i) - goodx(i)) < 1e-3)
    }
    
    val qpSolverUb = new QpSolver(n, 0, false, None, None, true, true)
    val ub = Array.fill[Double](n)(0.25)
    qpSolverUb.updateUb(ub)
    
    val (statusUb, ubResult) = qpSolverUb.run(ata, atb.mul(-1).data)
    for(i <- 0 until n) {
      println(ubResult(i))
    }
    
    val equalityBuilder = new CSCMatrix.Builder[Double](1,n)
    for(i <- 0 until n) equalityBuilder.add(0, i, 1)
    val qpSolverEq = new QpSolver(n, 0, false, Some(equalityBuilder.result), None, true, true)
    qpSolverEq.updateUb(Array.fill[Double](n)(1.0))
    qpSolverEq.updateEquality(Array[Double](0.5))
    val (statusEq, eqResult) = qpSolverEq.run(ata, atb.mul(-1).data)
    var sum: Double = 0.0
    for(i <- 0 until n) {
      println(eqResult(i))
      sum = sum + eqResult(i)
    }
    println("Sum " + sum + " status " + statusEq)
    assert(abs(sum - 0.5) < 1e-4)
    
    /* Test L1 constraints with OWLQN */
    val owlqn = new OWLQN[DenseVector[Double]](10, 4, 1.0)
    def optimizeWithOWLQN(init: DenseVector[Double]) : Boolean = {
      val f = new DiffFunction[DenseVector[Double]] {
        def calculate(x: DenseVector[Double]) = {
          ((math.pow(norm(x - 3.0,2),2)),(x * 2.0) - 6.0)
        }
      }
      val result = owlqn.minimize(f,init)
      for(i <- 0 until result.length) println(result(i))
      
      val closeish = norm(result - 2.5,2) < 1E-4
      closeish
    }
    assert(optimizeWithOWLQN(DenseVector.fill(10, 0.0)))
    
    /* Generate the Qp with L1 constraint */
    println("Generate Qp with L1 constraint")
    val l1 = 10
    val inequalityBuilder = new CSCMatrix.Builder[Double](2*l1, 2*l1)
    for (i <- 0 until l1) {
      inequalityBuilder.add(2*i, i, -1)
      inequalityBuilder.add(2*i, l1 + i, -1)
      inequalityBuilder.add(2*i + 1, i, 1)
      inequalityBuilder.add(2*i + 1, l1 + i, -1)
    }
    val qpSolverL1 = new QpSolver(l1, l1, true, None, Some(inequalityBuilder.result))
    val Hl1 = Array.fill[Double](l1)(2.0)
    val fl1 = Array.fill[Double](2*l1)(-6.0)
    for(i <- 0 until l1) fl1.update(l1 + i, 1.0)
    val (statusL1, resultL1) = qpSolverL1.run(Hl1, fl1)
    for(i <- 0 until l1) println(resultL1(i))
  }
}