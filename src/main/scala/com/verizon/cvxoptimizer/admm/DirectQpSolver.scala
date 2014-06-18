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
class DirectQpSolver(H: DoubleMatrix, alpha: Double, rho: Double,
  A: Option[Array[Double]] = None) {
  val n = H.rows

  def updateGram(row: Int, col: Int, value: Double) {
    if (row < 0 || row >= n)
      throw new IllegalArgumentException("ProximalQpSolver updateHessian row is out of range")
    if (col < 0 || col >= n)
      throw new IllegalArgumentException("ProximalQpSolver updateHessian col is out of range")
    H.data(row * n + col) = value
  }

  val MAX_ITER = 1000
  val ABSTOL = 1e-8
  val RELTOL = 1e-4
  
  for(i <- 0 until H.rows) {
    val diag = H.get(i, i)
    H.put(i, i, diag + rho)
  }
  
  val R = Decompose.cholesky(H)
  val Rtrans = R.transpose()
  
  var z = DoubleMatrix.zeros(n, 1)
  var u = DoubleMatrix.zeros(n, 1)
  
  var zOld = DoubleMatrix.zeros(n, 1)
  var xHat = DoubleMatrix.zeros(n, 1)
  var scale = DoubleMatrix.zeros(n, 1)

  var residual = DoubleMatrix.zeros(n, 1)
  var s = DoubleMatrix.zeros(n, 1)

  //Default is same as dposv: One cholesky factorization followed by forward-backward solves
  var proximal = 0

  def setProximal(p: Int): DirectQpSolver = {
    proximal = p
    this
  }

  def solve(q: DoubleMatrix, lb: Option[DoubleMatrix] = None, ub: Option[DoubleMatrix] = None): DoubleMatrix = {
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
      zOld.fill(0).addi(z).muli(1-alpha)
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
        case 5 => Proximal.proxL1(z.data, rho)
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
          
      if (residualNorm < epsPrimal && sNorm < epsDual) {
    	  println("DirectQpSolver converged in iterations " + k)
    	  return x
      }
      k = k + 1
    }
    println("DirectQpSolver MAX ITER reached convergence failure call ECOS")
    x
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
    val qpSolver = new DirectQpSolver(H,alpha,rho)
    val lb = DoubleMatrix.zeros(problemSize,1)
    val ub = DoubleMatrix.zeros(problemSize,1).addi(10)
    
    val qpStart = System.currentTimeMillis()
    val result = qpSolver.solve(f)
    val qpTime = System.currentTimeMillis() - qpStart
    
    assert(result.subi(3.0).norm2() < 1E-4)
    
    println("dim " + problemSize + " bfgs " + bfgsTime + " dposv " + dposvTime + " directqp " + qpTime)
  }
}