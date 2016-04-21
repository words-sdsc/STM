package optimizeLBFGS 

import breeze.linalg._ //{DenseVector, DenseMatrix, sum}
import breeze.numerics.{exp, log, sqrt}
import breeze.optimize.{LBFGS, DiffFunction}

object likelihood { 
  
  def lhoodFunction(beta: DenseMatrix[Double], doc_ct: DenseVector[Double], mu: DenseVector[Double], 
      siginv: DenseMatrix[Double]): DiffFunction[DenseVector[Double]] = {
    
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
        val part1 = betas*( doc_ct / sum(betas(::,*)).toDenseVector) - (expeta *(sum(doc_ct)/sum(expeta)))
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
 
//end likelihood
}