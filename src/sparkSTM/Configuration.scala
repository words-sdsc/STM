package sparkSTM

class Configuration {
  
  //dim
  var dim$K : Int = 0
  var dim$V : Int = 0
  var dim$A : Int = 0
  var dim$N : Int = 0

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
   
}