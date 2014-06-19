package com.verizon.cvxoptimizer.admm

import com.verizon.cvxoptimizer.proximal.Proximal
import org.jblas.DoubleMatrix
import org.jblas.Decompose
import org.jblas.Solve
import scala.math.min
import scala.math.max
import scala.math.sqrt
import breeze.optimize.LBFGS
import breeze.linalg.DenseVector
import breeze.optimize.DiffFunction
import breeze.linalg.norm

/*
 * Available with BSD license due to use of Proximal functions from BSD licensed code
 * 
 * Author @ Debasish Das
 * 
 * Reference: http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
 * 
 * Extension for proximal operators that show up in Matrix Factorization using Alternating Least Squares
 *  
 */

/*
 * Proximal operator and ADMM based Primal-Dual QP Solver 
 * 
 * It solves problem that has the following structure
 * 
 * 1/2 x*'Hx + f*'x
 * s.t 	x >= 0
 * 	 lb <= x <= ub
 *   ax = b
 *   ax <= b
 *   
 * When constraints are much less than variables it should behave better than ECOS 
 * We have to see what's the threshold when larger constraints will break the primal solver	
 */
class DirectQpSolver(n: Int,
    lb: Option[DoubleMatrix] = None, ub: Option[DoubleMatrix] = None,
    A: Option[Array[Double]] = None, alpha: Double = 1.0, rho: Double = 1.0) {
	
  val wsH = DoubleMatrix.zeros(n, n)
  
  val MAX_ITER = 1000
  val ABSTOL = 1e-8
  val RELTOL = 1e-4

  var z = DoubleMatrix.zeros(n, 1)
  var u = DoubleMatrix.zeros(n, 1)

  var zOld = DoubleMatrix.zeros(n, 1)
  var xHat = DoubleMatrix.zeros(n, 1)
  var scale = DoubleMatrix.zeros(n, 1)

  var residual = DoubleMatrix.zeros(n, 1)
  var s = DoubleMatrix.zeros(n, 1)
  
  var lambda = 1.0
  var proximal = 0
  
  def setLambda(lambda: Double) : DirectQpSolver = {
    this.lambda = lambda
    this
  }
  
  //TO DO : This can take a proximal function as input
  def setProximal(prox: Int) : DirectQpSolver = {
    this.proximal = prox
    this
  }
  
  def updateGram(row: Int, col: Int, value: Double) {
    if (row < 0 || row >= n) {
      throw new IllegalArgumentException("DirectQpSolver row out of bounds for gram matrix update")
    }
    if (col < 0 || col >= n) {
      throw new IllegalArgumentException("DirectQpSolver column out of bounds for gram matrix update")
    }
    wsH.put(row, col, value)
  }
  //Default is same as dposv: One cholesky factorization followed by forward-backward solves
  
  def solve(q: DoubleMatrix): DoubleMatrix = {
    //Dense cholesky factorization
    val R = Decompose.cholesky(wsH)
    val Rtrans = R.transpose()
    
    z.fill(0)
    u.fill(0)
    
    residual.fill(0)
    s.fill(0)
    
    var k = 0
    
    //Memory for x and tempR are allocated by Solve.solve calls
    //TO DO : See how this is implemented in breeze, why a workspace can't be used
    var x: DoubleMatrix = null
    var scaledR: DoubleMatrix = null
    
    while (k < MAX_ITER) {
      //scale = rho*(z - u) - q
      scale.fill(0)
      scale.addi(z).subi(u).muli(rho).subi(q)

      // x = R \ (R' \ scale)
      //Step 1 : scale * y = R'
      scaledR = Solve.solve(Rtrans, scale)

      //Step 2 : y * x = R
      x = Solve.solve(R, scaledR)

      //z-update with relaxation

      //zold = (1-alpha)*z
      //x_hat = alpha*x + zold
      zOld.fill(0).addi(z).muli(1 - alpha)
      xHat.fill(0).addi(x).muli(alpha).addi(zOld)

      //zold = z
      zOld.fill(0).addi(z)

      //z = xHat + u
      z.fill(0).addi(xHat).addi(u)

      //Apply proximal operator

      //Pick the correct proximal operator based on options
      //We will test the following

      //0. no proximal
      //1. proxBound
      //2. proxPos
      //3. prox + Linear
      //4. proxPos + Linear
      //5. proxL1
      proximal match {
        case 0 =>
        case 1 =>
          if (lb == None && ub == None)
            throw new IllegalArgumentException("DirectQpSolver proximal operator on box needs lower and upper bounds")
          Proximal.projectBox(z.data, lb.get.data, ub.get.data)
        case 2 => Proximal.projectPos(z.data)
        case 3 => if (A != None) Proximal.proxLinear(z.data, rho, A.get)
        case 4 => if (A != None) Proximal.proxLp(z.data, rho, A.get)
        case 5 => Proximal.proxL1(z.data, lambda*rho)
      }

      //z has proximal(x_hat)

      //Dual (u) update
      u.addi(xHat.subi(z))

      //Convergence checks
      //history.r_norm(k)  = norm(x - z);
      residual.fill(0).addi(x).subi(z)
      val residualNorm = residual.norm2()

      //history.s_norm(k)  = norm(-rho*(z - zold));
      s.fill(0).addi(z).subi(zOld).muli(-rho)
      val sNorm = s.norm2()

      //TO DO : Make sure z.muli(-1) is actually needed in norm calculation
      //residual = -z
      residual.fill(0).addi(z).muli(-1)
      //s = rho*u
      s.fill(0).addi(u).muli(rho)

      val epsPrimal = sqrt(n) * ABSTOL + RELTOL * max(x.norm2(), residual.norm2())
      val epsDual = sqrt(n) * ABSTOL + RELTOL * s.norm2()

      if (residualNorm < epsPrimal && sNorm < epsDual) return x
      k = k + 1
    }
    println("DirectQpSolver MAX ITER reached convergence failure call ECOS")
    x
  }
  
  def solve(H: DoubleMatrix, q: DoubleMatrix): DoubleMatrix = {
    for (i <- 0 until H.rows)
      for(j <- 0 until H.columns) {
        val h = H.get(i,j)
        if(i==j) wsH.put(i, j, rho + h)
        else wsH.put(i, j, h)
      }
    solve(q)
  }
}

object DirectQpSolver {
  def main(args: Array[String]) {
    if (args.length < 1) {
      println("Usage: QpSolver 100")
      println("Test QpSolver with a simple quadratic function of dimension 100")
      sys.exit(1)
    }

    val problemSize = args(0).toInt
    println("Test DirectQpSolver, breeze LBFGS and jblas solvePositive with problemSize " + problemSize)

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
    val jblasf = DoubleMatrix.zeros(problemSize, 1).add(-6.0)

    val dposvStart = System.currentTimeMillis()
    val dposvResult = Solve.solvePositive(jblasH, jblasf).data
    val dposvTime = System.currentTimeMillis() - dposvStart

    val H = DoubleMatrix.eye(problemSize).mul(2.0)
    val f = DoubleMatrix.zeros(problemSize, 1).add(-6)
    val alpha = 1.0
    val rho = 1.0

    val qpSolver = new DirectQpSolver(problemSize)
    val qpStart = System.currentTimeMillis()
    val result = qpSolver.solve(H, f)
    val qpTime = System.currentTimeMillis() - qpStart

    assert(result.subi(3.0).norm2() < 1E-4)

    println("dim " + problemSize + " bfgs " + bfgsTime + " dposv " + dposvTime + " directqp " + qpTime)
    
    val n = 5
    val ata = new DoubleMatrix(Array(
      Array(4.377, -3.531, -1.306, -0.139, 3.418),
      Array(-3.531, 4.344, 0.934, 0.305, -2.140),
      Array(-1.306, 0.934, 2.644, -0.203, -0.170),
      Array(-0.139, 0.305, -0.203, 5.883, 1.428),
      Array(3.418, -2.140, -0.170, 1.428, 4.684)))

    val atb = new DoubleMatrix(Array(-1.632, 2.115, 1.094, -1.025, -0.636))

    val goodx = Array(0.13025, 0.54506, 0.2874, 0.0, 0.028628)

    val qpSolverPos = new DirectQpSolver(n).setProximal(2)
    val posResult = qpSolverPos.solve(ata, atb.muli(-1))

    for (i <- 0 until n) {
      println(posResult.get(i) + " " + goodx(i))
      assert(Math.abs(posResult.get(i) - goodx(i)) < 1e-4)
    }
    
    val goodBounds: DoubleMatrix = DoubleMatrix.zeros(n, 1)
    goodBounds.put(0,0,0.0)
    goodBounds.put(1,0,0.25000000000236045)
    goodBounds.put(2,0,0.2499999999945758)
    goodBounds.put(3,0,0.0)
    goodBounds.put(4,0,0.0)

    val lb = DoubleMatrix.zeros(problemSize, 1)
    val ub = DoubleMatrix.zeros(problemSize, 1).addi(0.25)
    val qpSolverBounds = new DirectQpSolver(n,Some(lb), Some(ub)).setProximal(1)
    val boundsResult = qpSolverBounds.solve(ata, atb)
    println("Bounds result check " + (boundsResult.subi(goodBounds).norm2() < 1e-4))
    
    val qpSolverL1 = new DirectQpSolver(problemSize).setProximal(5)
    val l1Results = qpSolverL1.solve(H, f)
    println("L1 result check " + (l1Results.subi(2.5).norm2() < 1e-3))
    
    //Lp formulation
    //x1 + x2 + .. + xr <= 0, x1 >= 0, x2>= 0, ... , xr >= 0
    //Still have to generate a golden for this test
    /*
    val A = Array.fill[Double](problemSize)(1.0)
    val qpSolverLp = new DirectQpSolver(problemSize, None, None, Some(A)).setProximal(4)
    val lpResults = qpSolverLp.solve(H, f)
    println("Lp results " + lpResults)
    */
    
    println("Spark tests")
    
    val Htest = DoubleMatrix.zeros(1,1)

    Htest.put(0, 0, 1.123621)
    val ftest = DoubleMatrix.zeros(1,1)
    ftest.put(0, 0, -0.521311)
    
    val testSolver = new DirectQpSolver(1)
    val testResult = testSolver.solve(Htest, ftest)
    println(testResult)
    println(Solve.solvePositive(Htest, ftest.mul(-1)))
    
    val testResult1 = testSolver.solve(Htest, ftest)
    println(testResult1)
  }
}