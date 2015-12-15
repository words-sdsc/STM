package optimizeLBFGS

object likelihood {
  
  def main(args: Array[String]) 
  {
  import breeze.linalg._
  import breeze.numerics.{exp, log}
  import breeze.optimize._

  println("Maximization of Likelihood Function")
      
  def lhoodFunction(beta: DenseMatrix[Double], doc_ct: DenseVector[Double], mu: DenseVector[Double], siginv: DenseMatrix[Double]): DiffFunction[DenseVector[Double]] = {
    
    // f = -(objective function)
    val f = new DiffFunction[DenseVector[Double]] {
      
      override def valueAt(eta: DenseVector[Double]): Double = {
        val expetaRow = new DenseMatrix(rows=1, cols=eta.length + 1, data=(eta.data ++ Array(0.0)).map(x => exp(x)))
        val ndoc   = sum(doc_ct)
        val part1  = ((log(expetaRow * beta))*doc_ct - ndoc*log(sum(expetaRow)))
        val diff   = eta - mu
        val part2  = (diff.t * siginv * diff) * 0.5
        part2 - part1(0)
        
      }
      
      override def gradientAt(eta: DenseVector[Double]): DenseVector[Double] = {
        val expeta = new DenseVector((eta.data ++ Array(0.0)).map(x => exp(x)))
        val betas = beta(::, *) :* expeta                                     //each column of beta .* expeta (elementwise)
        val part1 = beta*( doc_ct / sum(beta(::,*)).toDenseVector) - (expeta *(sum(doc_ct)/sum(expeta)))
        //explicit transpose of row-sum removed, toDenseVector transposes it
        val part2 = siginv*(eta - mu) 
        val part1s = new DenseVector(part1.data.dropRight(1))   //drop the last row
        part2 - part1s
      }
      
      def calculate(eta: DenseVector[Double]): (Double, DenseVector[Double]) = {
        val J = valueAt(eta)
        val grad = gradientAt(eta)
        (J, grad)
        }
     
    } //end f
    
    f
  }// end lhoodFunction
    
    val lbfgs = new LBFGS[DenseVector[Double]](maxIter = 10000)    
    
    
    //**************************** Test Fixture : Start *************************
    
    val ntopics = 20
    val vacobSize = 10
    
    val siginvp    = DenseMatrix.rand(ntopics-1,ntopics-1)
    val initialEta = DenseVector.rand(ntopics-1)
    val mup        = DenseVector.rand(ntopics-1)
    val betap      = DenseMatrix.rand(ntopics, vacobSize)
    val doc_ctp    = linspace(1.0, vacobSize, vacobSize)
        
    // for each document, get f and minimize it {  
            println(initialEta)
            val likeliHoodF  = lhoodFunction(betap, doc_ctp, mup, siginvp)    
            val newEta       = lbfgs.minimize(likeliHoodF, initialEta)
            println(newEta)
    // }    

    //**************************** Test Fixture : End *************************

  //end main
  } 
 
//end likelihood
}