package optimizeLBFGS

import org.json.simple.parser.JSONParser
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;
import java.io.FileReader
import java.util.HashMap

import breeze.linalg.{DenseVector, DenseMatrix, sum}
import breeze.numerics.{exp, log, sqrt}
import breeze.optimize.{LBFGS}

object tests {
  
  //**************************** Spark HDFS [Test] *************************
  
  def hessPhiBound(x: String) = 
  {
    x match {
          case "hessian" => test_hessian()
          case "phi" => test_phi()
          case "bound" => test_bound()
    } 
    
          def test_hessian() = {
            print("\n||||||||: Testing Hessian")
          }
          
          def test_phi() = {
            print("\n||||||||: Testing Phi")
          }
          
          def test_bound() = {
            print("\n||||||||: Testing Bound")
          }
  }
  
  def sparkhdfs() = {
            import org.apache.spark.{SparkContext, SparkConf}
            import scala.math.random
      
            val conf = new SparkConf().setAppName("Spark Pi").setMaster("local")
            val spark = new SparkContext(conf)
            val textFile = spark.textFile("hdfs://localhost:9000/whereIsHadoop.txt")
            val counts = textFile.flatMap(line => line.split(" "))
                 .map(word => (word, 1))
                 .reduceByKey(_ + _)
            counts.saveAsTextFile("hdfs://localhost:9000/CountResults_whereIsHadoop")
            spark.stop()
  }
  
  //**************************** Likelihood [Test] *************************
  def llihood() = {
    println("Maximization of Likelihood Function") 
    
    //JSON import section * start
    try {
      val parser : JSONParser = new JSONParser();
      val reader : FileReader = new FileReader("/Users/aloksingh/Downloads/benchmark.JSON");
      val jsonObject :JSONObject = (parser.parse(reader)).asInstanceOf[JSONObject]
      
      val outp = jsonObject.get("benchmark").asInstanceOf[HashMap[String, JSONObject]]
      println(outp.keySet())
      
      var uuu = outp.get("hpbcpp").asInstanceOf[JSONObject].get("eta").asInstanceOf[JSONObject].get("lambda").asInstanceOf[JSONArray]
      var ggg = outp.get("grad").asInstanceOf[JSONArray] 
      
      
      //--
      var jsonoptiEtaSum = 1.0
      
      val jsonoptiEta : DenseVector[Double] = DenseVector.zeros[Double](uuu.size()) 
      val jsonoptiEtaNorm : DenseVector[Double] = DenseVector.zeros[Double](uuu.size()+1)     

       for(i <- 0 until uuu.size()) {
        jsonoptiEta(i) = (uuu.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
        jsonoptiEtaSum = jsonoptiEtaSum + exp(jsonoptiEta(i))
       }
      
      for(i <- 0 until uuu.size()) {
        jsonoptiEtaNorm(i) = exp(jsonoptiEta(i)) / jsonoptiEtaSum
      }
      jsonoptiEtaNorm(74) = 1.0 / jsonoptiEtaSum
      println("&&&&&&&&&&&&&&&&&&&&& opti eta json normalized &&&&&&&&&&&&&&&&&&&&&")
      println(jsonoptiEtaNorm)
      println("&&&&&&&&&&&&&&&&&&&&& opti eta json normalized &&&&&&&&&&&&&&&&&&&&&")
      //--
      
      val gradCheck : DenseVector[Double] = DenseVector.zeros[Double](ggg.size()) 
      val jsongradNorm : DenseVector[Double] = DenseVector.zeros[Double](ggg.size()+1)     

      var gradSum = 1.0
      for(i <- 0 until ggg.size()) {
        gradCheck(i) = (ggg.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
        gradSum = gradSum + exp(gradCheck(i))
       }
 
      
      for(i <- 0 until ggg.size()) {
        jsongradNorm(i) = exp(gradCheck(i)) / gradSum
      }
      jsongradNorm(74) = 1.0 / gradSum
      
      //read inputs from JSON file
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
      
      val likeliHoodF  = likelihood.lhoodFunction(betap, doc_ctp, mup, siginvp)
      val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=10) // maxIter = 1000000)
      val newEta       = lbfgs.minimize(likeliHoodF, initialEta)
      
      /*println("-----initialEta data output-------")
      println(initialEta)
      println(likeliHoodF.valueAt(initialEta))
      println(likeliHoodF.gradientAt(initialEta))
      
      println("-----newEta data output-------")
      println(newEta)
      println(likeliHoodF.valueAt(newEta))
      println(likeliHoodF.gradientAt(newEta))*/
      
      val newGradient = likeliHoodF.gradientAt(newEta)
       
      //normalizing newGradient @ newEta
      var newSum = 1.0
      for(i <- 0 until newGradient.length) {
        newSum = newSum + exp(newGradient(i))
      }
      
      val newgradNorm : DenseVector[Double] = DenseVector.zeros[Double](newGradient.length+1)     
      for(i <- 0 until newGradient.length) {
        newgradNorm(i) = exp(newGradient(i)) / newSum
      }
      newgradNorm(74) = 1.00 / newSum
      
      
      //normalizing newEta
      newSum = 1.0
      for(i <- 0 until newEta.length) {
        newSum = newSum + exp(newEta(i))
      }
      
      val newEtaNorm : DenseVector[Double] = DenseVector.zeros[Double](newEta.length+1)     
      for(i <- 0 until newEta.length) {
        newEtaNorm(i) = exp(newEta(i)) / newSum
      }
      newEtaNorm(74) = 1.00 / newSum
      
      
      println("normalized gradients*******************")
      println("Gradient parameter at new optimal eta from Scala")
      //println(newGradient)
      println("newEta normalized start")
      println(newEtaNorm)
      println("newEta normalized end")

      println(newgradNorm)
      println(newgradNorm.length)
      
      println("Gradient parameters at json output")
      //println(gradCheck)
      println(jsongradNorm.length)
      println(jsongradNorm)
      
      println("\n\ndifference b/w gradient parameter*******************")
      println(jsongradNorm - newgradNorm)
      println(sum(jsongradNorm - newgradNorm))

      
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
    
    
  }
}