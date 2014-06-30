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
 * X(i) is unbounded below, Set UB(i) = Inf if X(i) is unbounded above.
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
  
  var bounds = 0
  if (lbFlag) bounds += 1
  if (ubFlag) bounds += 1
  
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
  	Gquad = [fhalf', -1/sqrt(2),
         	-W, zerocolumn,
         	-fhalf', +1/sqrt(2)],
	hquad = [1/sqrt(2), zerocolumn, 1/sqrt(2)],
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
    if (nHessian != W.rows) {
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

  def solve(H: DoubleMatrix, f: Array[Double]): (Int, Array[Double]) = {
    if (diagonal) {
      throw new IllegalArgumentException("Qpsolver: digonal flag must be false for dense solve")
    }
    updateHessian(H)
    
    updateLinearObjective(f)
    
    val status = NativeECOS.solveSocp(c, G.rows, G.cols, G.data, G.colPtrs, G.rowIndices, hBuilder,
      Aeq.rows, Aeq.cols, Aeq.data, Aeq.colPtrs, Aeq.rowIndices, beqBuilder,
      linear, cones, x)
    
    (status, x.slice(0, nHessian))
  }

  def solve(Hdiag: Array[Double], f: Array[Double]): (Int, Array[Double]) = {
    if (!diagonal) {
      throw new IllegalArgumentException("QpSolver: diagonal flag must be true for sparse solve")
    }
    updateDiagonal(Hdiag)
    
    updateLinearObjective(f)
    
    val status = NativeECOS.solveSocp(c, G.rows, G.cols, G.data, G.colPtrs, G.rowIndices, hBuilder,
      Aeq.rows, Aeq.cols, Aeq.data, Aeq.colPtrs, Aeq.rowIndices, beqBuilder,
      linear, cones, x)
    
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
    val (statusDiag, qpResultDiag) = qpSolverDiag.solve(Hdiag, f)
    val qpDiagTime = System.currentTimeMillis() - qpDiagStart

    val qpSolver = new QpSolver(problemSize)
    val qpStart = System.currentTimeMillis()
    val (status, qpResult) = qpSolver.solve(H, f)
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
    
    val (statusBounds, qpResultBounds) = qpSolverBounds.solve(ata, atb.mul(-1).data)
    
    for (i <- 0 until n) {
      println(qpResultBounds(i) + " " + goodx(i))
      assert(Math.abs(qpResultBounds(i) - goodx(i)) < 1e-3)
    }
    
    val qpSolverUb = new QpSolver(n, 0, false, None, None, true, true)
    val ub = Array.fill[Double](n)(0.25)
    qpSolverUb.updateUb(ub)
    
    val (statusUb, ubResult) = qpSolverUb.solve(ata, atb.mul(-1).data)
    for(i <- 0 until n) {
      println(ubResult(i))
    }
    
    val equalityBuilder = new CSCMatrix.Builder[Double](1,n)
    for(i <- 0 until n) equalityBuilder.add(0, i, 1)
    val qpSolverEq = new QpSolver(n, 0, false, Some(equalityBuilder.result), None, true, true)
    qpSolverEq.updateUb(Array.fill[Double](n)(1.0))
    qpSolverEq.updateEquality(Array[Double](0.5))
    val (statusEq, eqResult) = qpSolverEq.solve(ata, atb.mul(-1).data)
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
    val (statusL1, resultL1) = qpSolverL1.solve(Hl1, fl1)
    for(i <- 0 until l1) println(resultL1(i))
    
    //TO DO : Add a testcase for L1 with dense H matrix
    
    println("Movielens test")
    println("Unbounded")
    val Hml = new DoubleMatrix(25, 25, 1244.661211, 885.841321, 899.989489, 932.895648, 877.170721, 836.475018, 904.400733, 853.370981, 1047.154258, 911.472203, 1004.418889, 1052.684672, 904.683742, 898.752730, 846.257240, 880.965938, 987.095303, 835.290194, 926.045883, 995.277621, 771.721543, 772.050185, 916.108833, 919.052137, 876.729461, 885.841321, 1145.670414, 848.745231, 839.949474, 815.990674, 782.962886, 771.258623, 791.912830, 939.233711, 849.336432, 893.230169, 942.649858, 848.374582, 817.650013, 767.540650, 798.900682, 939.461082, 776.694551, 813.977613, 906.668580, 731.382514, 718.822256, 858.858369, 865.998896, 793.789007, 899.989489, 848.745231, 1109.198294, 870.307368, 804.054911, 798.313080, 812.904963, 798.196046, 987.184288, 848.813247, 918.747185, 972.379634, 836.475997, 830.848025, 776.732930, 798.468257, 949.197865, 771.954054, 844.514021, 934.549742, 714.747793, 709.434575, 832.883741, 861.176778, 802.702918, 932.895648, 839.949474, 870.307368, 1185.873323, 850.495989, 828.199278, 815.778369, 855.709535, 1004.482635, 877.694364, 956.401361, 1001.956922, 911.555021, 880.032359, 802.514414, 869.037271, 904.178588, 819.007354, 882.023816, 973.570360, 781.771511, 717.951188, 920.487746, 875.004093, 812.477347, 877.170721, 815.990674, 804.054911, 850.495989, 1082.889103, 764.706225, 783.258270, 798.918797, 934.254442, 823.991178, 886.936654, 934.251085, 845.301478, 835.195949, 761.336034, 801.814763, 871.054422, 766.226641, 799.969291, 882.853531, 722.803309, 673.654165, 843.564472, 835.637433, 767.084530, 836.475018, 782.962886, 798.313080, 828.199278, 764.706225, 1030.453254, 751.279906, 784.487632, 901.911486, 790.508380, 867.174997, 908.247164, 826.876639, 816.707117, 736.108714, 795.472315, 859.816738, 727.746597, 811.499970, 863.042064, 697.923266, 636.115910, 811.396083, 795.887145, 770.207224, 904.400733, 771.258623, 812.904963, 815.778369, 783.258270, 751.279906, 1122.732950, 735.387458, 965.457072, 800.434884, 913.688664, 957.669829, 774.765297, 831.922470, 746.071727, 762.266589, 953.332558, 701.004806, 868.493004, 903.291715, 639.128663, 692.906338, 768.952840, 856.039554, 769.726777, 853.370981, 791.912830, 798.196046, 855.709535, 798.918797, 784.487632, 735.387458, 1067.459254, 895.938246, 801.498811, 891.407186, 905.904177, 865.254792, 832.543729, 760.606581, 841.057989, 840.274200, 778.227561, 817.761243, 879.868766, 748.216946, 656.781887, 850.486236, 814.386838, 755.196553, 1047.154258, 939.233711, 987.184288, 1004.482635, 934.254442, 901.911486, 965.457072, 895.938246, 1400.824752, 973.841940, 1077.711501, 1122.957392, 973.845338, 970.305097, 872.614898, 941.345253, 1081.502722, 871.281621, 998.299947, 1091.644096, 807.545882, 823.466420, 985.788248, 1001.753617, 923.116975, 911.472203, 849.336432, 848.813247, 877.694364, 823.991178, 790.508380, 800.434884, 801.498811, 973.841940, 1127.154700, 909.557045, 967.533540, 860.681911, 844.729459, 787.433410, 832.293838, 895.853526, 786.140374, 834.837372, 934.965270, 747.602639, 714.456792, 889.751871, 864.204947, 790.767159, 1004.418889, 893.230169, 918.747185, 956.401361, 886.936654, 867.174997, 913.688664, 891.407186, 1077.711501, 909.557045, 1286.917155, 1067.323971, 912.496825, 934.588007, 842.598541, 886.720330, 1022.609980, 835.479164, 948.511906, 1044.882401, 771.608958, 769.186134, 913.633712, 943.607940, 893.438031, 1052.684672, 942.649858, 972.379634, 1001.956922, 934.251085, 908.247164, 957.669829, 905.904177, 1122.957392, 967.533540, 1067.323971, 1392.008140, 965.176749, 978.416984, 892.246290, 959.345378, 1064.239900, 890.931306, 1000.607525, 1080.378329, 820.003292, 824.932331, 967.568892, 981.963499, 915.408427, 904.683742, 848.374582, 836.475997, 911.555021, 845.301478, 826.876639, 774.765297, 865.254792, 973.845338, 860.681911, 912.496825, 965.176749, 1173.179538, 869.884155, 808.593198, 870.599845, 881.842105, 814.142895, 855.410833, 950.013413, 808.456490, 703.446605, 934.727388, 855.966225, 812.922193, 898.752730, 817.650013, 830.848025, 880.032359, 835.195949, 816.707117, 831.922470, 832.543729, 970.305097, 844.729459, 934.588007, 978.416984, 869.884155, 1154.776492, 786.103933, 846.615903, 909.335274, 797.761454, 879.867447, 940.137270, 734.620394, 701.388486, 867.123126, 859.975145, 805.539361, 846.257240, 767.540650, 776.732930, 802.514414, 761.336034, 736.108714, 746.071727, 760.606581, 872.614898, 787.433410, 842.598541, 892.246290, 808.593198, 786.103933, 1002.484476, 768.388224, 832.989359, 737.081456, 779.535735, 854.898899, 691.478292, 644.686578, 801.010765, 789.647739, 733.117918, 880.965938, 798.900682, 798.468257, 869.037271, 801.814763, 795.472315, 762.266589, 841.057989, 941.345253, 832.293838, 886.720330, 959.345378, 870.599845, 846.615903, 768.388224, 1107.565113, 871.508396, 776.592009, 854.294568, 909.860416, 754.880472, 674.515444, 877.422803, 827.081749, 790.727727, 987.095303, 939.461082, 949.197865, 904.178588, 871.054422, 859.816738, 953.332558, 840.274200, 1081.502722, 895.853526, 1022.609980, 1064.239900, 881.842105, 909.335274, 832.989359, 871.508396, 1320.191790, 816.832622, 945.933513, 1016.373457, 737.516061, 786.473500, 885.188757, 944.327702, 874.515683, 835.290194, 776.694551, 771.954054, 819.007354, 766.226641, 727.746597, 701.004806, 778.227561, 871.281621, 786.140374, 835.479164, 890.931306, 814.142895, 797.761454, 737.081456, 776.592009, 816.832622, 1003.734452, 791.362150, 853.609073, 697.912412, 652.189201, 795.034436, 774.364882, 721.352180, 926.045883, 813.977613, 844.514021, 882.023816, 799.969291, 811.499970, 868.493004, 817.761243, 998.299947, 834.837372, 948.511906, 1000.607525, 855.410833, 879.867447, 779.535735, 854.294568, 945.933513, 791.362150, 1176.651156, 951.301184, 700.484670, 700.469987, 860.079680, 867.866514, 818.508476, 995.277621, 906.668580, 934.549742, 973.570360, 882.853531, 863.042064, 903.291715, 879.868766, 1091.644096, 934.965270, 1044.882401, 1080.378329, 950.013413, 940.137270, 854.898899, 909.860416, 1016.373457, 853.609073, 951.301184, 1314.702058, 800.700083, 771.237107, 967.597201, 945.620028, 876.117099, 771.721543, 731.382514, 714.747793, 781.771511, 722.803309, 697.923266, 639.128663, 748.216946, 807.545882, 747.602639, 771.608958, 820.003292, 808.456490, 734.620394, 691.478292, 754.880472, 737.516061, 697.912412, 700.484670, 800.700083, 943.185740, 579.991557, 814.942174, 727.857371, 683.349528, 772.050185, 718.822256, 709.434575, 717.951188, 673.654165, 636.115910, 692.906338, 656.781887, 823.466420, 714.456792, 769.186134, 824.932331, 703.446605, 701.388486, 644.686578, 674.515444, 786.473500, 652.189201, 700.469987, 771.237107, 579.991557, 856.381404, 689.395007, 721.830853, 647.458990, 916.108833, 858.858369, 832.883741, 920.487746, 843.564472, 811.396083, 768.952840, 850.486236, 985.788248, 889.751871, 913.633712, 967.568892, 934.727388, 867.123126, 801.010765, 877.422803, 885.188757, 795.034436, 860.079680, 967.597201, 814.942174, 689.395007, 1235.703834, 874.160742, 815.874871, 919.052137, 865.998896, 861.176778, 875.004093, 835.637433, 795.887145, 856.039554, 814.386838, 1001.753617, 864.204947, 943.607940, 981.963499, 855.966225, 859.975145, 789.647739, 827.081749, 944.327702, 774.364882, 867.866514, 945.620028, 727.857371, 721.830853, 874.160742, 1134.787259, 809.487849, 876.729461, 793.789007, 802.702918, 812.477347, 767.084530, 770.207224, 769.726777, 755.196553, 923.116975, 790.767159, 893.438031, 915.408427, 812.922193, 805.539361, 733.117918, 790.727727, 874.515683, 721.352180, 818.508476, 876.117099, 683.349528, 647.458990, 815.874871, 809.487849, 1034.201541)
    val fml = Array(-4878.593152, -4386.234229, -4460.392085, -4756.731215, -4354.481336, -4204.889415, -4416.351925, -4290.478906, -5312.345905, -4583.562028, -5021.710510, -5235.008167, -4601.387756, -4563.743787, -4163.355832, -4457.898625, -4843.623681, -4159.399627, -4647.099715, -5097.547640, -3889.942052, -3734.875516, -4773.758835, -4605.413425, -4296.966969)
    
    val ml = 25
    val qpML = new QpSolver(ml)
    val (statusMl, resultMl) = qpML.solve(Hml, fml)
    for(i <- 0 until ml) println(resultMl(i))
    
    println("L1")
  }
}