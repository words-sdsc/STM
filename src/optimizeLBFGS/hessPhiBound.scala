package optimizeLBFGS

object hessPhiBound 
{
  import breeze.linalg._
  import breeze.numerics.{exp, log, sqrt, abs}
  import breeze.optimize._
  
  // colvec eta, mat beta, colvec v doc_ct, colvec mu, mat siginv, colvec sigmaentropy
  def hessPhiBoundFunction(eta: DenseVector[Double], beta: DenseMatrix[Double], doc_ct: DenseVector[Double], 
      mu: DenseVector[Double], siginv: DenseMatrix[Double], sigmaentropy: DenseVector[Double]): 
      Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = {
   
      val expetaCol = new DenseVector(eta.data ++ Array(0.0)).map(x => exp(x))
      val theta     = expetaCol/sum(expetaCol)
      
      val EBa  = beta.copy
      val EBb  = EBa(::, *) :* expetaCol
      val EB   = EBb(*, ::) :* (sqrt(doc_ct) / sum(EBb(::, *)).toDenseVector)
      
      val hess = EB*EB.t - (theta*theta.t)*sum(doc_ct)

      diag(hess) -= sum(EB(*,::)) - theta * sum(doc_ct)
   
      val hessianM = hess(0 to hess.rows - 2, 0 to hess.cols - 2) + siginv
      
      var nu = DenseMatrix.rand[Double](hessianM.rows, hessianM.cols)
      
      //perform Cholesky
      try  { 
             nu = cholesky(hessianM).t    //upper triangular
           } catch {
      
             case e: Exception => { 
                         //update diagonal and recalculate Cholesky
                         val dvec = diag(hessianM)
                         
                         val pos = abs(hessianM)
                         
                         val magnitudes = sum(pos(*,::)) - abs(dvec)
                         
                         val dvecMin = DenseVector.rand[Double](dvec.length)
                         
                         for( i <- 0 until dvec.length ) { 
                               dvecMin(i) = min(dvec(i), magnitudes(i))
                         }
                         
                         diag(hessianM) := dvecMin 
                         nu = cholesky(hessianM).t
                          
                         }//end exception
             }//end catch
           
             val detTerm : Double = sum(log(diag(nu))) * -1.0
             val nuInverted       = inv(upperTriangular(nu))   //check do we need upper triangular
             val nuFinal          = nuInverted * nuInverted.t
             val diff             = eta - mu
             val bound            = log(theta.toDenseMatrix * beta)*doc_ct + detTerm - diff.t * siginv * diff * 0.5 - sigmaentropy
             //toDenseMatrix performs transpose so no explicit transpose needed
              
    (EB, (eta, nuFinal), bound(0))
    // phi = EB
  }
}