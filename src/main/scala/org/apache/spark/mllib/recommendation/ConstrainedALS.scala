/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * ALS matrix factorization with constraints formulated as a Quadratic Program and solved using a SOCP solver
*/
package org.apache.spark.mllib.recommendation

import org.apache.spark.rdd.RDD
import scala.util.Random
import scala.util.hashing.byteswap32
import scala.collection.mutable.ArrayBuffer
import org.apache.spark.Partitioner
import scala.collection.mutable.BitSet
import scala.Array.canBuildFrom
import org.apache.spark.SparkContext._
import org.apache.spark.util.Utils
import org.apache.spark.{Logging, HashPartitioner, Partitioner}
import org.apache.spark.storage.StorageLevel
import scala.math.{abs, sqrt}

class ConstrainedALS (
  private var numBlocks: Int,
  private var rank: Int,
  private var iterations: Int,
  private var lambda: Double,
  private var implicitPrefs: Boolean,
  private var alpha: Double,
  private var seed: Long = System.nanoTime()) extends ALS {

  /* Given an RDD of (partition(u), (u, p, r)) it computes userInLinks and userOutLinks */
  /* makeLinkRDDs takes RDD[blockId, (u, p, r)] and runs partitionBy on it using HashPartitioner on numBlocks 
     * 
     * val grouped = ratings.partitionBy(new HashPartitioner(numBlocks)) will group all blockIds together
     * 
     * All (u, p, r) that belongs to blockId_i will go to worker_i 
     * 
     * val links = grouped.mapPartitionsWithIndex will generate a RDD by applying a function to each partition of this RDD
     * while tracking the index of the original partition 
     * 
     * mapPartitionWithIndex will give a tuple of blockId and Array[(u, p, r)] associated with this blockId 
     * 
     */
  def createPartitioner(numBlocks: Int) = {
    new Partitioner {
      val numPartitions = numBlocks

      def getPartition(x: Any): Int = {       
        Utils.nonNegativeMod(byteswap32(x.asInstanceOf[Int]), numPartitions)
      }
    }
  }

  def makeInLinkBlock(numBlocks: Int, ratings: Array[Rating], partitioner: Partitioner): InLinkBlock = {
    val userIds = ratings.map(r => r.user).distinct.sorted
    val numUsers = userIds.length
    //userIdToPos HashMap
    val userIdToPos = userIds.zipWithIndex.toMap

    //generate ratings by product block
    //say we are running with 10 partitions
    //blockRatings[0]....blockRatings[9]
    val blockRatings = Array.fill(numBlocks)(new ArrayBuffer[Rating])
    
    //val blockRatings = Array.fill(numBlocks)(new ArrayBuffer[Rating])
    //Get each element from Array[Rating] and generate all products associated with each block 
    for (r <- ratings) {
      blockRatings(partitioner.getPartition(r.product)) += r
    }

    /*
       * In-link information for a user block. This includes the original user IDs of the elments within
       * this block, also an array of indices and ratings that specify which user in the block will be rated
       * by which products from each product block
       * 
       * users InLinkBlock, ratingsForBlock(b)(i) i-th product in productBlock b will have two arrays indices and ratings
       * indices are the index of users (by their index in this block) and ratings are the rating
       */
    val ratingsForBlock = new Array[Array[(Array[Int], Array[Double])]](numBlocks)

    for (productBlock <- 0 until numBlocks) {
      //Create array of (product, Seq(Rating)) ratings
      //blockRatings(productBlock) has all the Rating for this particular product Block
      //groupBy groups them by a Map[productId, ArrayBuffer[Rating]]
      //toArray will generate Array[(productId, ArrayBuffer[Rating])]

      val groupedRatings = blockRatings(productBlock).groupBy(_.product).toArray
      //sorted the groupedRatings by productId

      //groupedRating has sorted productIds

      //ratingsForBlock for each productBlock
      //Now we map the array to a new Array[Array[userIndices: Int], Array[rating: Double]]
      //productIds are sorted but the productId is not used in creating ratingsForBlock(productBlock)
      ratingsForBlock(productBlock) = groupedRatings.map {
        case (p, rs) =>
          (rs.view.map(r => userIdToPos(r.user)).toArray, rs.view.map(_.rating).toArray)
      }
    }
    //InLinkBlock(userIds, ratingsForBlock)
    ???
  }

  /*
     * Make the out-links table for a block of the users dataset give the list of (user, product, rating) values
     * for the users in that block
     */
  def makeOutLinkBlock(numBlocks: Int, ratings: Array[Rating], partitioner: Partitioner): OutLinkBlock = {
    val userIds = ratings.map(_.user).distinct.sorted
    val numUsers = userIds.length
    val userIdToPos = userIds.zipWithIndex.toMap
    val shouldSend = Array.fill(numUsers)(new BitSet(numBlocks))
    for (r <- ratings) {
      shouldSend(userIdToPos(r.user))(partitioner.getPartition(r.product)) = true
    }
    OutLinkBlock(userIds, shouldSend)
  }

  private def makeLinkRDDs(numBlocks: Int, ratings: RDD[(Int, Rating)], partitioner: Partitioner): (RDD[(Int, InLinkBlock)], RDD[(Int, OutLinkBlock)]) =
    {
      val grouped = ratings.partitionBy(new HashPartitioner(numBlocks))

      val links = grouped.mapPartitionsWithIndex((blockId, elements) => {
        val ratings = elements.map { _._2 }.toArray
        val inLinkBlock = makeInLinkBlock(numBlocks, ratings, partitioner)
        val outLinkBlock = makeOutLinkBlock(numBlocks, ratings, partitioner)
        Iterator.single((blockId, (inLinkBlock, outLinkBlock)))
      }, true)
      val inLinks = links.mapValues(_._1)
      val outLinks = links.mapValues(_._2)
      inLinks.persist(StorageLevel.MEMORY_AND_DISK)
      outLinks.persist(StorageLevel.MEMORY_AND_DISK)
      (inLinks, outLinks)
    }
  
  private def randomFactor(rank: Int, rand: Random): Array[Double] = {
    // Choose a unit vector uniformly at random from the unit sphere, but from the
    // "first quadrant" where all elements are nonnegative. This can be done by choosing
    // elements distributed as Normal(0,1) and taking the absolute value, and then normalizing.
    // This appears to create factorizations that have a slightly better reconstruction
    // (<1%) compared picking elements uniformly at random in [0,1].
    val factor = Array.fill(rank)(abs(rand.nextGaussian()))
    val norm = sqrt(factor.map(x => x * x).sum)
    factor.map(x => x / norm)
  }
  
  override def run(ratings: RDD[Rating]): MatrixFactorizationModel = {
    val sc = ratings.context

    val numBlocks = if (this.numBlocks == -1) {
      math.max(sc.defaultParallelism, ratings.partitions.size / 2)
    } else {
      this.numBlocks
    }

    val partitioner = createPartitioner(numBlocks)

    /* partition(u), (u, p, r) */
    val ratingsByUserBlock = ratings.map { rating =>
      (partitioner.getPartition(rating.user), rating)
    }

    /* partition(p), (p, u, r) */
    val ratingsByProductBlock = ratings.map { rating =>
      (partitioner.getPartition(rating.product),
        Rating(rating.product, rating.user, rating.rating))
    }

    //generate userInLinks, userOutLinks
    val (userInLinks, userOutLinks) = makeLinkRDDs(numBlocks, ratingsByUserBlock, partitioner)

    //shuffle rating.user, rating.product to rating.product, rating.user so that same makeLinkRDDs function can be used
    val (productInLinks, productOutLinks) = makeLinkRDDs(numBlocks, ratingsByProductBlock, partitioner)

    // Initialize user and product factors randomly, but use a deterministic seed for each
    // partition so that fault recovery works
    val seedGen = new Random(seed)
    val seed1 = seedGen.nextInt()
    val seed2 = seedGen.nextInt()

    /* generate randomFactor for each user of size rank */
    var users = userOutLinks.mapPartitionsWithIndex { (index, itr) => /* index is blockId, iter is Iterator[(userId, OutLinkBlock)] */
      val rand = new Random(byteswap32(seed1 ^ index))
      itr.map {
        case (x, y) => (x, y.elementIds.map(_ => randomFactor(rank, rand)))
      }
    }
    
    /* generate randomFactor for each product of size rank */
    var products = productOutLinks.mapPartitionsWithIndex { (index, itr) =>
      val rand = new Random(byteswap32(seed2 ^ index))
      itr.map {
        case (x, y) =>
          (x, y.elementIds.map(_ => randomFactor(rank, rand)))
      }
    }

    for (iter <- 1 to iterations) {
      //perform ALS update
      logInfo("Re-computing I given U (Iteration %d/%d)".format(iter, iterations))
      products = updateFeatures(users, userOutLinks, productInLinks, partitioner, rank, lambda, alpha)
      logInfo("Re-computing U given I (Iteration %d/%d)".format(iter, iterations))
      users = updateFeatures(products, productOutLinks, userInLinks, partitioner, rank, lambda, alpha)
    }

    ???
  }

  /*
   * Compute the user feature vectors given the current products (or vice-versa).
   * 
   * products: RDD[(productId, features)] features is Array[Array[Double]] Array[Double] is of rank r
   * OutLinkBlock computes all the productIds in each block, the productIds are stored in OutLinkBlock.elementIds
   * 
   * For this productId, 
   */
  def updateFeatures(
    products: RDD[(Int, Array[Array[Double]])],
    productOutLinks: RDD[(Int, OutLinkBlock)],
    userInLinks: RDD[(Int, InLinkBlock)],
    partitioner: Partitioner,
    rank: Int,
    lambda: Double,
    alpha: Double): RDD[(Int, Array[Array[Double]])] = {
    val numBlocks = products.partitions.size
    //We are updating the features for user, to generate user features we need product features...for each productId we ha
    productOutLinks.join(products)
    ???
  }
}