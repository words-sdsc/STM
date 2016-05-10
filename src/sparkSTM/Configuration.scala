package sparkSTM
import breeze.linalg.DenseMatrix

class Configuration {
  
  //dim
  var dim$K : Int = 0
  var dim$V : Int = 0
  var dim$A : Int = 0
  var dim$N : Int = 0
  var dim$wcounts$x : List[Int] = null
  
  //covariates
  var covariates$betaindex = null

  //kappa
  var kappa$LDAbeta: scala.Any    = null
  var kappa$Interactions: scala.Any = null
  
  //convergence
  var convergence$threshold: Double = 0.0
  var convergence$max_iterations: Int = 0
  
  //init
  var init$verbose : Boolean = false
  var init$mode = null
  var init$nits = null
  var init$alpha = null
  var init$eta = null
  var init$burnin = null
  var init$s = null
  var init$p = null
  var init$d_group_size = null
  
  //covariance
  var covariance : DenseMatrix[Double] = null  //used by update_mu
  var sigprior : Double = 0.0 //used by update_sigma
}