package sparkSTM

import breeze.linalg.{DenseMatrix}

object start {
  
  def main(args: Array[String]) 
  { 
    
    //tests.test_hessPhiBound()
    //tests.test_sparkhdfs()
    //tests.test_likelihood() 
 
    val model = new STMModel()
    val settings = new Configuration()
    val documents : List[DenseMatrix[Double]] = null
    model.initialize(documents, settings)
  
  }
}