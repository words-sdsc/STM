package optimizeLBFGS

import java.io.FileReader
import java.util.HashMap
 
import org.json.simple.parser.JSONParser
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;


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
    
    val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=11)    //maxIter = 10000 removed tolerance = 1E-28
    
    
    //**************************** Test Fixture : Start *************************
    
    //JSON import section * start
    try {
      val parser : JSONParser = new JSONParser();
      val reader : FileReader = new FileReader("/Users/aloksingh/Downloads/benchmark.JSON");
      val jsonObject :JSONObject = (parser.parse(reader)).asInstanceOf[JSONObject]
      
      val outp = jsonObject.get("benchmark").asInstanceOf[HashMap[String, JSONObject]]
      println(outp.keySet())
      
      var ggg = outp.get("grad").asInstanceOf[JSONArray] 
      val gradCheck : DenseVector[Double] = DenseVector.zeros[Double](ggg.size())     
      for(i <- 0 until ggg.size()) {
        gradCheck(i) = (ggg.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
       }
      
      val f  = jsonObject.get("inputs").asInstanceOf[HashMap[String, JSONObject]]
      println(f.keySet())
       
      var temp = f.get("mu").asInstanceOf[JSONArray]     
      val mup : DenseVector[Double] = DenseVector.zeros[Double](temp.size())     
      for(i <- 0 until temp.size()) {
        mup(i) = temp.get(i).toString().toDouble
      }
      
      var temp2 = f.get("doc_ct").asInstanceOf[JSONArray]     
      val doc_ctp : DenseVector[Double] = DenseVector.zeros[Double](temp2.size())     
      for(i <- 0 until temp2.size()) {
        doc_ctp(i) = temp2.get(i).toString().toDouble
      }
      
      var temp3 = f.get("eta").asInstanceOf[JSONArray]     
      val initialEta : DenseVector[Double] = DenseVector.zeros[Double](temp3.size())     
      for(i <- 0 until temp3.size()) {
        initialEta(i) = temp3.get(i).toString().toDouble
      }
      
      var temp4 = f.get("beta").asInstanceOf[JSONArray]       
      val betap : DenseMatrix[Double] = DenseMatrix.zeros[Double](75,108) 
      for(i <- 0 until temp4.size()) {
            var ccc = temp4.get(i).asInstanceOf[JSONArray]     
            for(j <- 0 until ccc.size()) {
              betap(i,j) = ccc.get(j).toString().toDouble
            }
      }
            
      var temp5 = f.get("siginv").asInstanceOf[JSONArray]       
      val siginvp : DenseMatrix[Double] = DenseMatrix.zeros[Double](74,74)     
      for(i <- 0 until temp5.size()) {
            var eee = temp5.get(i).asInstanceOf[JSONArray]     
            for(j <- 0 until eee.size()) {
              siginvp(i,j) = eee.get(j).toString().toDouble
            }
      }
            
      /*val initialEta = DenseVector.rand(ntopics-1)
      val mup        = DenseVector.rand(ntopics-1)
      val betap      = DenseMatrix.rand(ntopics, vacobSize)
      val doc_ctp    = linspace(1.0, vacobSize, vacobSize) */
      
      val likeliHoodF  = lhoodFunction(betap, doc_ctp, mup, siginvp)    
      val newEta       = lbfgs.minimize(likeliHoodF, initialEta)
      
      println("-----initialEta data output-------")
      println(initialEta)
      println(likeliHoodF.valueAt(initialEta))
      println(likeliHoodF.gradientAt(initialEta))
      
      println("-----newEta data output-------")
      println(newEta)
      println(likeliHoodF.valueAt(newEta))
      println(likeliHoodF.gradientAt(newEta))
      
      /*println("-----test data output-------")
      println(gradCheck)  
      println("-----diff of grad output-------")
      println(likeliHoodF.gradientAt(newEta) - gradCheck)
      println(likeliHoodF.gradientAt(newEta) - likeliHoodF.gradientAt(newEta)) */
      
    } catch {
          case e: Exception => println(e.getMessage)
    }
    
    //JSON import section end
    
    /*val ntopics = 20
    val vacobSize = 10
    val siginvp    = DenseMatrix.rand(ntopics-1,ntopics-1)
    val initialEta = DenseVector.rand(ntopics-1)
    val mup        = DenseVector.rand(ntopics-1)
    val betap      = DenseMatrix.rand(ntopics, vacobSize)
    val doc_ctp    = linspace(1.0, vacobSize, vacobSize) */
     

    //**************************** Test Fixture : End *************************

  //end main
  } 
 
//end likelihood
}