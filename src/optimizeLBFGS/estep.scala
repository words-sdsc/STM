package optimizeLBFGS

object estep {
  import breeze.linalg._
  import breeze.optimize.LBFGS
  import breeze.numerics.{exp, log, sqrt}

  //import org.apache.spark.mllib.linalg


  // colvec eta, mat beta, colvec v doc_ct, colvec mu, mat siginv, colvec sigmaentropy
  def evaluate(documents: List[DenseVector[Double]], betaIndex: DenseVector[Int], 
      updateMu: Boolean, beta: List[DenseMatrix[Double]], lambdaOld: List[DenseVector[Double]],
      mu: DenseMatrix[Double], sigma: DenseMatrix[Double],
      verbose: Boolean): 
      Tuple4[DenseMatrix[Double], DenseVector[DenseMatrix[Double]],
        DenseVector[Double], DenseMatrix[Double] ] = {
    
    //declare the constants
    val V = beta(1).cols
    val K = beta(1).rows
    val N = documents.length
    val A = beta.length
    
    //initialization
    val sigma_g =  DenseMatrix.zeros[Double](K-1,K-1)
    val beta_g  =  DenseVector.zeros[DenseMatrix[Double]](A)
    for(i <- 0 until beta_g.size) { 
      beta_g(i) = DenseMatrix.zeros[Double](K,V) 
    }
    
    val bound  = DenseVector.zeros[Double](N)
    val lambda = DenseMatrix.zeros[Double](K-1, N) //every COL is a lambda of len K-1, transpose it later
    DenseVector.zeros[Double](N)
     
    //preprocessing of common components
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
    
    
    //iterate over documents 1...N
    for(j <- 0 until N) {
      val doc    = documents(j)
      val words  = documents(j).toArray.map(_.toInt).toIndexedSeq  // words vs doc
      val aspect = betaIndex(j)
      val init   = lambdaOld(j)
      //if(updateMu)  //check if we need this
      val mu_j = mu(::,j) 
      val beta_j = beta(aspect)(::, words).toDenseMatrix
      
      //infer the document
      val results : Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = 
      
        logisticnormal(init, mu_j, siginv, beta_j, doc, sigmaentropy)
       
      //update the global values
      
      //sigma
      sigma_g += results._2._2
      //beta
      beta_g(aspect)(::, words) += results._1.asInstanceOf[Matrix[Double]]
      //bound
      bound(j) = results._3
      //lambda
      lambda(::, j) := results._2._1
      
      /*Dimensions:
       * val initialEta = DenseVector.rand(ntopics-1)
       * val mup        = DenseVector.rand(ntopics-1)
       * val betap      = DenseMatrix.rand(ntopics, vacobSize)
       * val doc_ctp    = linspace(1.0, vacobSize, vacobSize) */
       
    }
     
    //combine and return(list(sigma, beta, bound, lambda))
      
    (sigma_g, beta_g, bound, lambda)
  } //end of function evaluate
  
  def logisticnormal(eta: DenseVector[Double], mu: DenseVector[Double], 
        siginv: DenseMatrix[Double], beta: DenseMatrix[Double], doc_ct: DenseVector[Double], 
        sigmaentropy: Double): 
      Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = {
    
      val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=11)    //maxIter = 10000 removed tolerance = 1E-28
      val likeliHoodF  = likelihood.lhoodFunction(beta, doc_ct, mu, siginv)    
      val newEta       = lbfgs.minimize(likeliHoodF, eta)
      
      hessPhiBound.evaluate(newEta, beta, doc_ct, mu, siginv, sigmaentropy )
  }
  
}