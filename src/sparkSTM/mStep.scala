/**
 * Structural Topic Model by Molly Roberts, Brandon Stewart, Dustin Tingley
 * Website: http://structuraltopicmodel.com/
 * R Package: https://github.com/bstewart/stm
 * 
 * Scala:
 * @author  Alok Singh (a1singh@ucsd.edu)
 */

package sparkSTM

import breeze.linalg.{DenseMatrix, DenseVector, diag, sum, inv,*}
import breeze.numerics.{abs}

object mStep {

    def update_sigma(nu:DenseMatrix[Double],lambda:DenseMatrix[Double],
        mu:DenseMatrix[Double],sigprior: Double) :DenseMatrix[Double] = {
      val M = lambda - mu.t
      val covariance = M.t * M
      var sigma = (covariance+nu) / (1.0*lambda.rows)
      sigma = diag(diag(sigma))*sigprior + sigma * (1-sigprior)
      sigma
    }
     
    //***  used by update_mu() ***f() = Variational Linear Regression with a Half-Cauchy hyperprior 
    def vbreg(Y: DenseVector[Double], X:DenseMatrix[Double]) : DenseVector[Double] = {
      //b0=1, d0=1 taken out of function signature and replaced with constant
      val Xcorr  = X.t * X
      val XYcorr = X.t * Y
      val an = (1.0 + X.rows)/2.0
      val D  = X.cols
      val N  = X.rows
      val cn = X.cols
      var dn : Double = 1.0
      var Ea : Double = cn/dn
      var ba : Double = 1.0
      var w = DenseVector.zeros[Double](D)
      var precision : Double = 1.0
      var converge = 1000.0
      
      while(converge > 0.0001) {
        val wOld = w.copy
        val ppmat : DenseMatrix[Double] = diag(DenseVector.fill(D){Ea})
        ppmat(0,0) = 0.0
        val invV = ppmat + Xcorr :* precision
        val V : DenseMatrix[Double] = inv(invV)
        w = (V :* precision) * XYcorr
        
        val sse = sum((X * w - Y):^ 2.0)
        val bn = (sse + sum(diag(Xcorr * V))) * 0.5 + ba
        precision = an/bn
        ba = 1.0/(precision + 1.0) //b0=1
        
        val da = 2.0 / (Ea + 1.0)  //d0=1
        val diagV = diag(V)
        val w1 = w(1 to -1)
        dn = 2.0 * da + (w1.t * w1 + sum(diagV(1 to -1)))
        Ea = cn / dn
        converge = sum(abs(w - wOld))
      }
      
      w
    }
      
    def update_mu(lambda: DenseMatrix[Double], covar: DenseMatrix[Double]) : Tuple2[DenseMatrix[Double], DenseMatrix[Double]] = {
      //gamma = apply f(on every col of lambda) and cbind the result to get gamma matrix
      val gamma  : DenseMatrix[Double]  = lambda(::, * ).map( Y => vbreg(Y, covar) )
      val mu     : DenseMatrix[Double]  = (covar * gamma).t
      (mu, gamma)
    }
     
    def update_beta(betaList : DenseVector[DenseMatrix[Double]]) : DenseMatrix[Double] = {
       val B = betaList.foldLeft(DenseMatrix.zeros[Double](betaList(0).rows, betaList(0).cols))(_ + _)
       B(*, ::).map{ row => row / sum(row) } //row normalized
    }
}