package com.verizon.jecos

/*
 import org.scalatest.WrapWith
import org.specs.runner.ScalaTestSuite
import org.scalatest.junit._
import org.scalatest.prop._
*/

import breeze.optimize.LBFGS
import breeze.linalg.DenseVector
import breeze.optimize.DiffFunction
import breeze.linalg.norm
/*
import org.scalatest.FunSuite
@WrapWith(classOf[ScalaTestSuite])
class QpSolverTest extends FunSuite with Checkers {
  test("optimize a simple multivariate gaussian using QpSolver and compare it with Breeze LBFGS") {
    val lbfgs = new LBFGS[DenseVector[Double]](100, 4)
    def optimizeThis(init: DenseVector[Double]) = {
      val f = new DiffFunction[DenseVector[Double]] {
        def calculate(x: DenseVector[Double]) = {
          (norm((x-3.0):^ 2.0, 1), (x * 2.0) - 6.0)
        }
      }
      val result = lbfgs.minimize(f, init)
      norm(result - 3.0, 2) < 1E-10
    }
    val init = DenseVector.zeros[Double](100)
    init(0 until init.length by 2) := -1.2
    init(1 until init.length by 2) := 1.0
    assert(optimizeThis(init))
  }
}
*/