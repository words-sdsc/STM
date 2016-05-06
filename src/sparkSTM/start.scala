package sparkSTM

import breeze.linalg.{DenseMatrix}

object start {
  
  def main(args: Array[String]) 
  {  
    val model    = new STMModel()
    val settings = new Configuration()
    
    // ↳configure the settings object
     
    // ↳fill documents
    val documents : List[DenseMatrix[Double]] = null
    
    // ↳initialize the STM model
    model.initialize(documents, settings)
  
    val ntokens = settings.dim$wcounts$x.sum
    val betaindex = settings.covariates$betaindex
    
    // ↳run STM until convergence criteria or max iterations reached
    // model.runSTM(documents, betaIndex, updateMu, beta, lambdaOld, mu, sigma, verbose)
    
    // ↳construct output
    model.contruct_output()
    
    // ***other tests below this line ***
    //tests.test_hessPhiBound()
    //tests.test_sparkhdfs()
    //tests.test_likelihood() 
  }
}