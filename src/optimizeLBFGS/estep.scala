package optimizeLBFGS

object estep {
  import breeze.linalg._
  import breeze.optimize.LBFGS
  import breeze.numerics.{log}

  //import org.apache.spark.mllib.linalg


  def evaluate(documents: List[DenseMatrix[Double]], betaIndex: DenseVector[Int], 
      updateMu: Boolean, beta: List[DenseMatrix[Double]], lambdaOld: DenseMatrix[Double],
      mu: DenseMatrix[Double], sigma: DenseMatrix[Double],
      verbose: Boolean): 
      Tuple4[DenseMatrix[Double], DenseVector[DenseMatrix[Double]],
        DenseVector[Double], DenseMatrix[Double] ] = {
    
    //declare the constants
    val V = beta(1).cols
    val K = beta(1).rows
    val N = documents.length
    val A = beta.length //several docs could have same beta, hence A and N could be different
    
    //initialization
    val sigma_g =  DenseMatrix.zeros[Double](K,K)
    val beta_g  =  DenseVector.zeros[DenseMatrix[Double]](A)
    for(i <- 0 until beta_g.size) { 
      beta_g(i) = DenseMatrix.zeros[Double](K,V) 
    }
    
    val bound  = DenseVector.zeros[Double](N)
    val lambda = DenseMatrix.zeros[Double](K-1, N) //every COL is a lambda of len K-1, 
                                                   //transpose it before returning
     
     
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
      val wordIndices   = (documents(j)(0, ::).t).data.map(_.toInt).toIndexedSeq 
      //1st row of doc = indices of words
      val aspect = betaIndex(j)
      val init   = lambdaOld(j,::).t //matrix N x K-1, each row is a lambda for a doc
      //if(updateMu)  //check if we need this
      val mu_j   = mu(::,j) 
      val beta_j = beta(aspect)(::, wordIndices).toDenseMatrix
      
      //infer the document
      val results : Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = 
        logisticnormal(init, mu_j, siginv, beta_j, doc(1,::).t, sigmaentropy)
      //2nd row of doc = count of words

      //update the global values
      
      //sigma
      sigma_g += results._2._2
      //beta
      beta_g(aspect)(::, wordIndices) += results._1.asInstanceOf[Matrix[Double]]
      //bound
      bound(j) = results._3
      //lambda
      lambda(::, j) := results._2._1 //each COL is a lambda
       
    }
           
    (sigma_g, beta_g, bound, lambda.t)
    //transpose lambda so that each ROW is a lambda of length K-1
    
  } //end of function evaluate
  
  def logisticnormal(eta: DenseVector[Double], mu: DenseVector[Double], 
        siginv: DenseMatrix[Double], beta: DenseMatrix[Double], doc_ct: DenseVector[Double], 
        sigmaentropy: Double): 
      Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = {
    
      val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=11)    //maxIter = 10000 removed tolerance = 1E-28
      val likeliHoodF  = likelihood.lhoodFunction(beta, doc_ct, mu, siginv) 
                         //how does it get which doc_ct is for which word in dict ?
                         //--> adjust #rows in beta as per #rows in doc_ct (drop all zero rows) maintain order
      val newEta       = lbfgs.minimize(likeliHoodF, eta)
      
      hessPhiBound.evaluate(newEta, beta, doc_ct, mu, siginv, sigmaentropy )
  }
  
}