/**
 * Structural Topic Model by Molly Roberts, Brandon Stewart, Dustin Tingley
 * Website: http://structuraltopicmodel.com/
 * R Package: https://github.com/bstewart/stm
 * 
 * Scala:
 * @author  Alok Singh (a1singh@ucsd.edu)
 */

package sparkSTM

import breeze.linalg.{DenseMatrix}

object start {
  
  def main(args: Array[String]) 
  {  
    
    tests.test_initialization()
    
    
    
    
    /*
    val model  = new STMModel()
    
    val config = new Configuration()
    // ↳configure the settings object
     
    val documents : List[DenseMatrix[Double]] = null
    // ↳fill documents
    
    model.initialize(documents, config)
    // ↳initialize the STM model
    
    val ntokens = config.dim$wcounts$x.sum
    val betaindex = config.covariates$betaindex
    
    //model.runSTM(documents, betaIndex, updateMu, beta, lambdaOld, mu, sigma, verbose)
    // ↳run STM until convergence criteria or max iterations reached
    
    model.contruct_output() */
    // ↳construct output
    
    
    
    // ***other tests below this line ***
    //tests.test_hessPhiBound()
    //tests.test_sparkhdfs()
    //tests.test_likelihood() 
  }
}