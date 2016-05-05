package sparkSTM


import org.apache.spark.rdd.RDD
import org.la4j.matrix.sparse.{CRSMatrix => SparseMatrix}
import breeze.linalg.{Axis, DenseMatrix, DenseVector, Matrix, diag, inv, sum, det, cholesky, all, *}
import breeze.optimize.LBFGS
import org.apache.spark.{SparkContext, SparkConf}
import breeze.numerics.log
import scala.{Iterator=>ScalaIterator, Double}
import java.util.Iterator
import scala.collection.mutable.ArrayBuffer

 
class STMModel {
  //GG = set of global parameters
  var beta_g  : DenseMatrix[Double] = null
  var mu_g    : Tuple2[DenseMatrix[Double], DenseMatrix[Double]] = null
  var sigma_g : DenseMatrix[Double] = null
  
  
  /*[function]********************* initialize the STM model **********************/
  def initialize(documents: List[DenseMatrix[Double]], settings: Configuration) = {
     println("Initializing the STM model (Spectral mode) ...")
     
     var K : Int = settings.dim("K")
     val V : Int = settings.dim("V")
     val N : Int = settings.dim("N")
     val A : Int = settings.dim("A")
     val verbose : Boolean = settings.init$verbose
     
           if(K >= V) {
             println("Error: Topics >= Words in vocab")
             System.exit(1)
           }
     //**************************ıllıllı ıllıllı**************************
     //{•------» (1) Prep the Gram matrix «------•}
     //**************************ıllıllı ıllıllı**************************
     val mat : SparseMatrix = spectral.docsToSparseMatrix(documents)
     var wprob : DenseVector[Double] = spectral.colSums(mat) 
     wprob /= sum(wprob)
     var Q : DenseMatrix[Double] = spectral.gram(mat)
     var Qsums : DenseVector[Double] = sum(Q(*, ::)) //sum of each row
     var keep: Seq[Int] = null
     var whichzero: Seq[Int] = null
     
     if(!all(Qsums)) {
       //there are some zeros
       whichzero = spectral.whichZeros(Qsums) //which indices have zero in the input vector
       keep = (0 to Qsums.length-1 toList) diff whichzero.toList
       Q = Q.delete(whichzero, Axis._0)
       Q = Q.delete(whichzero, Axis._1)
       Qsums = spectral.dropelements(Qsums, whichzero)
       wprob = spectral.dropelements(wprob, whichzero)
     } 
     //divide every col by row sum
     Q = Q(::,*) :/ Qsums
     
     
     //**************************ıllıllı ıllıllı************************** 
     //{•------» (2) anchor words «------•}
     //**************************ıllıllı ıllıllı************************** 
     var anchor : List[Int] = null
     if(verbose) println("Finding anchor words...")
     if(K!=0) {
        anchor = spectral.fastAnchor(Q, K, verbose)
     } else {
        anchor = spectral.tsneAnchor(Q)
        K = anchor.length
     }
     
     
     //**************************ıllıllı ıllıllı**************************
     //{•------» (3) recoverKL «------•}
     //**************************ıllıllı ıllıllı**************************
     if(verbose) println("Recovering Initialization ")
     var beta0 = spectral.recoverL2(Q, anchor, wprob, verbose) //[?] $A
     
     if(keep != null) { 
       beta0 = spectral.refillZeros(K, V, beta0, keep, whichzero)
     }
     
     //**************************ıllıllı ıllıllı**************************
     //{•------» (4) generate other parameters «------•}
     //**************************ıllıllı ıllıllı**************************
     val mu = DenseMatrix.zeros[Double](K-1, 1)
     val sigma = diag(DenseVector.fill(K-1){20.0})
     val lambda = DenseMatrix.zeros[Double](N, K-1)
     if(verbose) println("Initialization complete")
     
     var beta : List[DenseMatrix[Double]] = List(beta0)
     for (i <- 1 to A-1) {
       beta = beta0 :: beta
     }
     
     val model = (mu, sigma, beta, lambda)
  }
  
  def kappa_init() = {
    
  }
  
  /*[function]********************* run STM over a set of documents **********************/
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
                ScalaIterator((sigma_pt, beta_pt, bound_pt, lambda_pt))
        
      }//end mapPartitions
    
    //PP = aggregation of results from each partition in 'metricsFromPartitions'
      
      val collect_sigma : DenseMatrix[Double] = partitionsMetrics.map(_._1).treeAggregate(DenseMatrix.zeros[Double](K,K))(_ += _, _ += _)
      val collect_lambda: DenseMatrix[Double] = DenseMatrix.vertcat(partitionsMetrics.map(_._4).flatMap(list => list).collect().map(_.toDenseMatrix): _*)
      val collect_bound : DenseVector[Double] = DenseVector(partitionsMetrics.map(_._3).reduce{(prev, item) => prev ++ item}.toArray) 
      val collect_beta  : DenseVector[DenseMatrix[Double]] = partitionsMetrics.map(_._2).reduce((vec1,vec2) => vec1 :+ vec2)
      
      //update global parameters GG using these aggregates PP
      
        val covar_dummy : DenseMatrix[Double] = null                   // [ ? covar ]
        this.mu_g     = mStep.update_mu(collect_lambda, covar_dummy)
        
        val sigprior_dummy : Double = 0.5                              // [ ? is this a double]
        this.sigma_g  = mStep.update_sigma(collect_sigma, collect_lambda, mu_g._1, sigprior_dummy)
        
        this.beta_g   = mStep.update_beta(collect_beta)
        
    //unpersist broadcasted variables
    betaBc.unpersist()
    
    this
    
  } //end runSTM
  
  
  /*[function]********************* Input:sigma, Output:siginverse, sigmaentropy **********************/
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
  
  
  /*[function]********************* infer single doc **********************/
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