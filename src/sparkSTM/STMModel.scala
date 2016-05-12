package sparkSTM

import org.apache.spark.rdd.RDD
import org.la4j.matrix.sparse.{CRSMatrix => SparseMatrix}
import breeze.linalg.{Axis, DenseMatrix, DenseVector, Matrix, diag, inv, sum, det, cholesky, all, *}
import breeze.optimize.LBFGS
import org.apache.spark.{SparkContext, SparkConf}
import breeze.numerics.{abs, log}
import scala.{Iterator=>ScalaIterator, Double}
import java.util.Iterator
//import scala.collection.mutable.ArrayBuffer

 
class STMModel {
  //variables used in testing
  var fastanchorL : List[Int] = null
  var recoverL2M  : DenseMatrix[Double] = null
  
  //-------------------------
  
  var beta_init : DenseVector[DenseMatrix[Double]] = null
  
  //GG = set of global parameters
  var beta_g  : DenseMatrix[Double]  = null
  var mu_g    : Tuple2[DenseMatrix[Double], DenseMatrix[Double]] = (null, null) //mu_g = (mu, gamma )
  var sigma_g : DenseMatrix[Double]  = null
  var lambda_g: DenseMatrix[Double]  = null
  
  var settings: Configuration        = null
  var modelConvergence   : Convergence = null //info whether this model has converged
  
  /*[function]********************* initialize the STM model **********************/
  def initialize(documents: List[DenseMatrix[Double]], config: Configuration) = {
     println("Initializing the STM model (Spectral mode) ...")
     
     this.settings = config
     var K : Int = settings.dim$K
     val V : Int = settings.dim$V
     val N : Int = settings.dim$N
     val A : Int = settings.dim$A
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
       //println("//there are some zeros")
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
        if(settings.testmode) {
          this.fastanchorL = anchor
          System.out.println("\nfastAnchor calculated ... updated model.fastanchorL")
        }
     } else {
        anchor = spectral.tsneAnchor(Q)
        K = anchor.length
     }
     
     
     //**************************ıllıllı ıllıllı**************************
     //{•------» (3) recoverKL «------•}
     //**************************ıllıllı ıllıllı**************************
     if(verbose) println("Recovering Initialization ")
     var beta0 = spectral.recoverL2(Q, anchor, wprob, verbose)
     
     
     if(keep != null) { 
       //if(verbose) println("keep != null refilling zeros")
       beta0 = spectral.refillZeros(K, V, beta0, keep, whichzero)
       if(settings.testmode) {
          this.recoverL2M = beta0
          System.out.println("recoverL2 calculated ...updated model.recoverL2M")
       }
     }
     
     //**************************ıllıllı ıllıllı**************************
     //{•------» (4) generate other parameters «------•}
     //**************************ıllıllı ıllıllı**************************
     val mu = DenseMatrix.zeros[Double](K-1, 1)
     val sigma = diag(DenseVector.fill(K-1){20.0})
     val lambda = DenseMatrix.zeros[Double](N, K-1)
     if(verbose) println("Initialization complete")
     
     //assign beta for each aspect
     val beta = DenseVector.fill(A){beta0}  //List of Matrices
     
     
     //model = (mu, sigma, beta, lambda) //mu_g = (mu, gamma)
     this.mu_g       =(mu, null)                //mu_g = (mu, gamma)
     this.sigma_g    =sigma
     this.beta_init  =beta
     this.lambda_g   =lambda
  } 

  
  /*[function]********************* run STM over a set of documents **********************/
  def runSTM(documents: List[DenseMatrix[Double]], betaIndex: DenseVector[Int], 
      updateMu: Boolean, beta: DenseVector[DenseMatrix[Double]], lambdaOld: DenseMatrix[Double],
      mu: DenseMatrix[Double], sigma: DenseMatrix[Double],
      verbose: Boolean) = { 
    
    //Spark
    val conf   = new SparkConf().setAppName("Spark Pi").setMaster("local")
    val spark  = new SparkContext(conf)
    //val numPartitions = 10 //number of partitions
    val documentsRDD: RDD[(DenseMatrix[Double], Int)] = spark.parallelize(documents.zipWithIndex)//, numPartitions)
      
    //get copy of global parameters
    
    val V = beta(1).cols
    val K = beta(1).rows
    val A = beta.length //several docs could have same beta, hence A and N could be different
    val (siginv, sigmaentropy) = valuesFromSigma(sigma: DenseMatrix[Double])
    
    //broadcast any large data structures to each node if they are needed by every doc
    val betaBc      = spark.broadcast(beta)
    
    // check : is broadcasting of below variables efficient
    val betaIndexBc = spark.broadcast(betaIndex)
    val lambdaOldBc = spark.broadcast(lambdaOld)
    val muBc        = spark.broadcast(mu)
    
    while(! modelConvergence.stopits ) {
      
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
                    val aspect  = betaIndexBc.value(indx)
                    val init_   = lambdaOldBc.value(indx,::).t //matrix N x K-1, each row is a lambda for a doc
                    //if(updateMu)  //check if we need this
                    val mu_     = muBc.value(::,indx) 
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
      
      this.mu_g     = mStep.update_mu(collect_lambda, settings.covariance)      
      this.sigma_g  = mStep.update_sigma(collect_sigma, collect_lambda, mu_g._1, settings.sigprior)
      this.beta_g   = mStep.update_beta(collect_beta)
    
      //check CONVERGENCE
      checkModelConvergence(collect_bound)
      
      if(!modelConvergence.stopits & verbose) {
        println("iteration :"+modelConvergence.its)
      }
      
    } // end of while loop
    
    //unpersist broadcasted variables
    betaBc.unpersist()
    betaIndexBc.unpersist()
    lambdaOldBc.unpersist()
    muBc.unpersist()
    
    if(verbose) println("all iterations finished...")
    //println("call construct_output()...") 
    
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
  
  /*[function]********************* checks if this model has converged ***********************
  	 *updates member variable "modelConvergence" of this object*/
  
  def checkModelConvergence(bound: DenseVector[Double]) = {
    //settings: Configuration and state: Convergence are local variables of this model
    
    val verbose = settings.init$verbose
    val emtol = settings.convergence$threshold
    val maxits = settings.convergence$max_iterations
    
    if(modelConvergence == null) { modelConvergence = new Convergence() }
    
    modelConvergence.bound ::= sum(bound)    //0th position holds latest value in the bound: List
    
    if(modelConvergence.its > 1) {
      val oldB = modelConvergence.bound(1)
      val newB = modelConvergence.bound(0)
      val check = (newB-oldB)/abs(oldB)
            if(check < emtol) {
                  modelConvergence.converged = true
                  modelConvergence.stopits = true
                  if(verbose) println("Model Converged")
                  //return
            } else if(modelConvergence.its == maxits) {
                  modelConvergence.stopits = true
                  if(verbose) println("Model Terminated Before Convergence Reached")
                  //return
            }
    } 
    
    modelConvergence.its += 1
    //nothing is returned; 'modelConvergence' is member variable of this class
  }
  
    
  def contruct_output() = {
    
  }
  
  def kappa_init() = {
    
  }
  
}