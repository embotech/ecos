package com.verizon.cvxoptimizer.admm

import com.verizon.cvxoptimizer.proximal.Proximal
import org.jblas.DoubleMatrix
import org.jblas.Decompose
import org.jblas.Solve
import scala.math.min
import scala.math.max
import scala.math.sqrt

/*
 * Available with BSD license due to use of Proximal functions from BSD licensed code
 * 
 * Author @ Debasish Das
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
  val ABSTOL = 1e-4
  val RELTOL = 1e-2

  val R = Decompose.cholesky(H.add(DoubleMatrix.eye(n).mul(rho)))
  val Rtrans = R.transpose()

  var z = DoubleMatrix.zeros(n, 1)
  var u = DoubleMatrix.zeros(n, 1)

  var zOld = DoubleMatrix.zeros(n, 1)
  var xHat = DoubleMatrix.zeros(n, 1)
  var scale = DoubleMatrix.zeros(n, 1)

  var residual = DoubleMatrix.zeros(n, 1)
  var s = DoubleMatrix.zeros(n, 1)

  var proximal = 1

  def setProximal(p: Int): DirectQpSolver = {
    proximal = p
    this
  }

  def solve(q: DoubleMatrix, lb: Option[DoubleMatrix], ub: Option[DoubleMatrix]): DoubleMatrix = {
    z.fill(0)
    u.fill(0)

    zOld.fill(0)
    xHat.fill(0)
    scale.fill(0)

    residual.fill(0)
    s.fill(0)

    var k = 0

    //Memory for x and tempR are allocated by Solve.solve calls
    //TO DO : See how this is implemented in breeze, why a workspace can't be used
    var x: DoubleMatrix = null
    var tempR: DoubleMatrix = null

    while (k < MAX_ITER) {
      //scale = rho*(z - u) - q
      scale.copy(z).subi(u).muli(rho).subi(q)
      // x = R \ (R' \ scale)
      //Step 1 : scale * y = R'
      tempR = Solve.solve(scale, Rtrans)
      //Step 2 : y * x = R
      x = Solve.solve(tempR, R)

      //z-update with relaxation

      //x_hat = alpha*x + (1-alpha)*zold
      zOld.copy(z).muli(1 - alpha)
      xHat.copy(x).muli(alpha).addi(zOld)

      //Apply proximal operator
      zOld.copy(z)

      //Pick the correct proximal operator based on options
      //We will test the following

      //1. no proximal
      //2. proxBound
      //3. proxPos
      //4. prox + Linear
      //5. proxPos + Linear
      //6. proxL1
      proximal match {
        case 2 =>
          if (lb == None || ub == None)
            throw new IllegalArgumentException("DirectQpSolver proximal operator on box needs lower and upper bounds")
          Proximal.projectBox(xHat.data, lb.get.data, ub.get.data)
        case 3 => Proximal.projectPos(xHat.data)
        case 4 => if (A != None) Proximal.proxLinear(xHat.data, rho, A.get)
        case 5 => if (A != None) Proximal.proxLp(xHat.data, rho, A.get)
        case 6 => Proximal.proxL1(xHat.data, rho)
      }
      z.copy(xHat)
      
      //Dual (u) update
      u.addi(xHat.subi(z))

      //Convergence checks
      //history.r_norm(k)  = norm(x - z);
      val residualNorm = residual.copy(x).subi(z).norm2()
      //history.s_norm(k)  = norm(-rho*(z - zold));
      val sNorm = s.copy(z).subi(zOld).muli(-rho).norm2()

      //TO DO : Make sure z.muli(-1) is actually needed in norm calculation
      residual.copy(z).muli(-1)
      s.copy(u).muli(rho)

      val epsPrimal = sqrt(n) * ABSTOL + RELTOL * max(x.norm2(), residual.norm2())
      val epsDual = sqrt(n) * ABSTOL + RELTOL * s.norm2()
      if (residualNorm < epsPrimal && sNorm < epsDual) x
      k = k + 1
    }
    println("DirectQpSolver MAX ITER reached without convergence, call Interior Point Solver")
    x
  }
}