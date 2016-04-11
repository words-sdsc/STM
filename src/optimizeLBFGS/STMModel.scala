package optimizeLBFGS

import org.apache.spark.rdd.RDD
import org.apache.spark.mllib.linalg.{Vector}
import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.optimize.LBFGS
import org.apache.spark.{SparkContext, SparkConf}
import breeze.linalg.{Matrix, diag, inv, sum, det, cholesky}
import breeze.numerics.{log}

// Update: - the 'estep' is implemented in this class
class STMModel {
  //GG = set of global parameters
  
    
  /* ********************* run STM over a set of documents ********************* */
  def runSTM(documents: List[DenseMatrix[Double]], betaIndex: DenseVector[Int], 
      updateMu: Boolean, beta: List[DenseMatrix[Double]], lambdaOld: DenseMatrix[Double],
      mu: DenseMatrix[Double], sigma: DenseMatrix[Double],
      verbose: Boolean) : STMModel = { 
    
    //Spark
    val conf   = new SparkConf().setAppName("Spark Pi").setMaster("local")
    val spark  = new SparkContext(conf)
    val documentsRDD: RDD[(DenseMatrix[Double], Int)] = spark.parallelize(documents.zipWithIndex)
      
    //get copy of global parameters
    

    
    //declare the constants
    val V = beta(1).cols
    val K = beta(1).rows
    val N = documents.length
    val A = beta.length //several docs could have same beta, hence A and N could be different
    val (siginv, sigmaentropy) = getSigma(sigma: DenseMatrix[Double])
    
    //broadcast any large data structures to each node if they are needed by every doc
    
    
    //partition the set of documents
    val metricsFromPartitions: RDD[(DenseMatrix[Double], DenseVector[DenseMatrix[Double]], 
                                    DenseVector[Double], DenseMatrix[Double])] =   
      documentsRDD.mapPartitions 
      {
        docs => 
                
                //TP = set of variables specific to this partition
                val sigma_pt   =  DenseMatrix.zeros[Double](K,K) 
                val beta_pt    =  DenseVector.zeros[DenseMatrix[Double]](A)  //List of Matrices
                for(i <- 0 until beta_pt.size) { beta_pt(i) = DenseMatrix.zeros[Double](K,V) } 
                val bound_pt   = DenseVector.zeros[Double](N)
                val lambda_pt  = DenseMatrix.zeros[Double](K-1, N)  
                //every COL is a lambda of len K-1, transpose it before returning 
                
                docs.foreach { 
                  case (doc: DenseMatrix[Double], indx: Int) =>  
                       
                    val wordIndices   = (doc(0, ::).t).data.map(_.toInt).toIndexedSeq //1st row of doc = indices of words
                    val aspect = betaIndex(indx)
                    val init   = lambdaOld(indx,::).t //matrix N x K-1, each row is a lambda for a doc
                    //if(updateMu)  //check if we need this
                    val mu_j   = mu(::,indx) 
                    val beta_j = beta(aspect)(::, wordIndices).toDenseMatrix
                    
                    /* infer this document via STM,
                       DD = set of results from STM for this single document */
                    val results = logisticnormal(init, mu_j, siginv, beta_j, doc(1,::).t, sigmaentropy)
                                                                              //doc(2nd row)=count of words
                    
                    // update partition accumulators TP with DD
                    beta_pt(aspect)(::, wordIndices) += results._1.asInstanceOf[Matrix[Double]] //two docs can have same aspect
                    lambda_pt(::, indx)     := results._2._1 //each COL is a lambda   
                    sigma_pt                += results._2._2
                    bound_pt(indx)           = results._3

                }

                // return iterator (each tuple is a set of accumulators from one partition)
                Iterator((sigma_pt, beta_pt, bound_pt, lambda_pt.t))
        
      }//end mapPartitions
    
    /*each tuple in 'metricsFromPartitions' is from one partition
      PP = aggregation of results from each partition in 'metricsFromPartitions'
      */
      
    
    //update global parameters GG using these aggregates PP
    
    
    //unpersist broadcasted variables
    
    
    this
    
  } //end runSTM
  
   
  /* *********************returns siginverse, sigmaentropy ********************* */
  def getSigma(sigma: DenseMatrix[Double]) : Tuple2[DenseMatrix[Double], Double] = {
      //pre-processing of common components
      var sigobj       = DenseMatrix.zeros[Double](1,1)
      var sigmaentropy = 0.0
      val siginv       = inv(sigma)

      try  { 
             sigobj       = cholesky(sigma).t        //upper triangular
             sigmaentropy = sum(log(diag(sigobj)))
           } catch {
             case e: Exception => { 
             sigmaentropy = 0.5 * (log(det(sigma).abs)) 
           } 
             
           }//end catch
    
      (siginv, sigmaentropy)
  }//end getSigma
  
  
  /* *********************infer single doc ********************* */
  def logisticnormal(eta: DenseVector[Double], mu: DenseVector[Double], 
        siginv: DenseMatrix[Double], beta: DenseMatrix[Double], doc_ct: DenseVector[Double], 
        sigmaentropy: Double): 
      Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = {
    
      val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=11)    //maxIter = 10000 removed tolerance = 1E-28
      val likeliHoodF  = likelihood.lhoodFunction(beta, doc_ct, mu, siginv) 
                         //how does it get which doc_ct is for which word in dict ?
                         //--> adjust #cols in beta = #rows in doc_ct (drop all zero rows& maintain order)
                         //--> beta*doc_ct
      val newEta       = lbfgs.minimize(likeliHoodF, eta)
      
      hessPhiBound.evaluate(newEta, beta, doc_ct, mu, siginv, sigmaentropy )
  }
}