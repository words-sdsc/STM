package optimizeLBFGS

object hessPhiBound {
  import breeze.linalg._
  import breeze.numerics.{exp, log, sqrt}
  import breeze.optimize._
  
  // col vec eta, mat beta, col v doc_ct, col v mu, mat siginv, col v sigmaentropy
  def lhoodFunction(eta: DenseVector[Double], beta: DenseMatrix[Double], doc_ct: DenseVector[Double], 
      mu: DenseVector[Double], siginv: DenseMatrix[Double], sigmaentropy: DenseVector[Double]): 
      Tuple3[DenseMatrix[Double],Tuple2[DenseVector[Double], DenseMatrix[Double]],Double] = {
   
      val expetaCol = new DenseVector(eta.data ++ Array(0.0)).map(x => exp(x))
      val theta     = expetaCol/sum(expetaCol)
      
      val EB1  = beta.copy
      val EB2 = EB1(::, *) :* expetaCol
      val EB = EB2(*, ::) :* (sqrt(doc_ct) / sum(EB2(::, *)).toDenseVector)
      
      val hess = EB*EB.t - (theta*theta.t)*sum(doc_ct)

      diag(hess) -= sum(EB(*,::)) - theta * sum(doc_ct)
   
      val hessian = hess(0 to hess.rows - 2, 0 to hess.cols - 2) + siginv
      
      val nu = DenseMatrix.rand[Double](hessian.rows, hessian.cols)
      
      
    (EB, (eta, nu), bound)
  }
}