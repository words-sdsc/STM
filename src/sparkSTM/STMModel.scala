package sparkSTM

import org.apache.spark.rdd.RDD
import breeze.linalg.{DenseMatrix, DenseVector, Matrix, diag, inv, sum, det, cholesky}
import breeze.optimize.LBFGS
import org.apache.spark.{SparkContext, SparkConf}
import breeze.numerics.log
import scala.Iterator

 
class STMModel {
  //GG = set of global parameters
  var beta_g  : DenseMatrix[Double] = null
  var mu_g    : Tuple2[DenseMatrix[Double], DenseMatrix[Double]] = null
  var sigma_g : DenseMatrix[Double] = null
  
  /* ********************* run STM over a set of documents ********************* */
  def runSTM(documents: List[DenseMatrix[Double]], betaIndex: DenseVector[Int], 
      updateMu: Boolean, beta: DenseVector[DenseMatrix[Double]], lambdaOld: DenseMatrix[Double],
      mu: DenseMatrix[Double], sigma: DenseMatrix[Double],
      verbose: Boolean) : STMModel = { 
    
    //Spark
    val conf   = new SparkConf().setAppName("Spark Pi").setMaster("local")
    val spark  = new SparkContext(conf)
    val numPartitions = 10 //number of partitions
    val documentsRDD: RDD[(DenseMatrix[Double], Int)] = spark.parallelize(documents.zipWithIndex, numPartitions)
      
    //get copy of global parameters
    
    val V = beta(1).cols
    val K = beta(1).rows
    val A = beta.length //several docs could have same beta, hence A and N could be different
    val (siginv, sigmaentropy) = valuesFromSigma(sigma: DenseMatrix[Double])
    
    //broadcast any large data structures to each node if they are needed by every doc
    val betaBc = spark.broadcast(beta)
    
    //partition the set of documents
    val partitionsMetrics: RDD[(DenseMatrix[Double], DenseVector[DenseMatrix[Double]], 
                                    List[Double], List[DenseVector[Double]] )] =   
      documentsRDD.mapPartitions 
      {
        docs => 
                //declare the constants
                val N = docs.length
                
                //TP = set of variables specific to this partition
                val sigma_pt   =  DenseMatrix.zeros[Double](K,K) 
                val beta_pt    =  DenseVector.zeros[DenseMatrix[Double]](A)  //List of Matrices
                for(i <- 0 until beta_pt.size) { beta_pt(i) = DenseMatrix.zeros[Double](K,V) } 
                var bound_pt   = List[Double]()                         //DenseVector.zeros[Double](N)
                var lambda_pt  = List[DenseVector[Double]]()            //DenseMatrix.zeros[Double](K-1, N)  
                //every COL is a lambda of len K-1, transpose it before returning 
                
                
                docs.foreach { 
                  case (doc: DenseMatrix[Double], indx: Int) =>  
                       
                    val wordIndices   = (doc(0, ::).t).data.map(_.toInt).toIndexedSeq //1st row of doc = indices of words
                    val aspect  = betaIndex(indx)
                    val init_   = lambdaOld(indx,::).t //matrix N x K-1, each row is a lambda for a doc
                    //if(updateMu)  //check if we need this
                    val mu_     = mu(::,indx) 
                    val beta_   = betaBc.value(aspect)(::, wordIndices).toDenseMatrix
                    
                    /* infer this document via STM,
                       DD = set of results from STM for this single document */
                    val results = logisticnormal(init_, mu_, siginv, beta_, doc(1,::).t, sigmaentropy)
                                                                              //doc(2nd row)=count of words
                    
                    // update partition accumulators TP with DD
                    beta_pt(aspect)(::, wordIndices) += results._1.asInstanceOf[Matrix[Double]] //two docs can have same aspect
                    lambda_pt = results._2._1 :: lambda_pt                                      //each COL is a lambda   
                    sigma_pt                += results._2._2
                    bound_pt = results._3 :: bound_pt
                }
               
                // return iterator (each tuple is a set of accumulators from one partition)
                Iterator((sigma_pt, beta_pt, bound_pt, lambda_pt))
        
      }//end mapPartitions
    
    //PP = aggregation of results from each partition in 'metricsFromPartitions'
      
      val collect_sigma : DenseMatrix[Double] = partitionsMetrics.map(_._1).treeAggregate(DenseMatrix.zeros[Double](K,K))(_ += _, _ += _)
      val collect_lambda: DenseMatrix[Double] = DenseMatrix.vertcat(partitionsMetrics.map(_._4).flatMap(list => list).collect().map(_.toDenseMatrix): _*)
      val collect_bound : DenseVector[Double] = DenseVector(partitionsMetrics.map(_._3).reduce{(prev, item) => prev ++ item}.toArray) 
      val collect_beta  : DenseVector[DenseMatrix[Double]] = partitionsMetrics.map(_._2).reduce((vec1,vec2) => vec1 :+ vec2)
      
      //update global parameters GG using these aggregates PP
      
        val covar_dummy : DenseMatrix[Double] = null                   // [ ? covar ]
        mu_g     = mStep.update_mu(collect_lambda, covar_dummy)
        
        val sigprior_dummy : Double = 0.5                              // [ ? is this a double]
        sigma_g  = mStep.update_sigma(collect_sigma, collect_lambda, mu_g._1, sigprior_dummy)
        
        beta_g   = mStep.update_beta(collect_beta)
        
    //unpersist broadcasted variables
    betaBc.unpersist()
    
    this
    
  } //end runSTM
  
  /* ******Input:sigma*********Output:siginverse, sigmaentropy ********************* */
  def valuesFromSigma(sigma: DenseMatrix[Double]) : Tuple2[DenseMatrix[Double], Double] = {
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
  }//end valuesFromSigma
  
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