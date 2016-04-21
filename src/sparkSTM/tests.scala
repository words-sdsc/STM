package sparkSTM

import org.json.simple.parser.JSONParser
import org.json.simple.JSONArray
import org.json.simple.JSONObject
import java.io.FileReader
import java.util.HashMap
import breeze.linalg.{DenseVector, DenseMatrix, sum}
import breeze.numerics.abs
import breeze.optimize.{LBFGS}

object tests {
  var jsonObject   : JSONObject    = null 
  var etaJSON, hessianJSON, phiJSON, boundJSON : JSONArray = null 
  
  //inputs       
  var mup : DenseVector[Double] = DenseVector.zeros[Double](74)   
  var doc_ctp : DenseVector[Double] = DenseVector.zeros[Double](108) 
  var initialEta : DenseVector[Double] = DenseVector.zeros[Double](74)   
  var betap : DenseMatrix[Double] = DenseMatrix.zeros[Double](75,108) 
  var siginvp : DenseMatrix[Double] = DenseMatrix.zeros[Double](74,74) 
  var sigmaentropyp: Double = 0.0
  
  //**************************** hessPhiBound [Test] *************************
  
  def test_hessPhiBound() = 
  {
    try {
      //read input from JSON file
      val parser : JSONParser = new JSONParser();
      val reader : FileReader = new FileReader("/Users/aloksingh/Downloads/benchmark2.JSON");
      jsonObject              = (parser.parse(reader)).asInstanceOf[JSONObject]  
      val outp = jsonObject.get("benchmark").asInstanceOf[HashMap[String, JSONObject]]
      //println(outp.keySet())
      
      fillInputs()
      
      //outputs from JSON for verification
      boundJSON = outp.get("hpbcpp").asInstanceOf[JSONObject].get("bound").asInstanceOf[JSONArray]
      etaJSON   = outp.get("hpbcpp").asInstanceOf[JSONObject].get("eta").asInstanceOf[JSONObject].get("lambda").asInstanceOf[JSONArray]
      hessianJSON = outp.get("hpbcpp").asInstanceOf[JSONObject].get("eta").asInstanceOf[JSONObject].get("nu").asInstanceOf[JSONArray]
      phiJSON = outp.get("hpbcpp").asInstanceOf[JSONObject].get("phis").asInstanceOf[JSONArray]
      
    } catch {
                  case e: Exception => println(e.getMessage)
            }
    
    test_all()
    
          def test_all() = {
            print("\n||||||||: Testing Bound") 
   
            //X : process using current implementation
            val likeliHoodF  = likelihood.lhoodFunction(betap, doc_ctp, mup, siginvp)
            val lbfgs = new LBFGS[DenseVector[Double]](tolerance = 1E-12, m=10) // maxIter = 1000000)
            val newEta       = lbfgs.minimize(likeliHoodF, initialEta) 
            
            val jsonoptiEta : DenseVector[Double] = DenseVector.zeros[Double](etaJSON.size())
            for(i <- 0 until etaJSON.size()) {
                    jsonoptiEta(i) = (etaJSON.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
            }
            
            val X = hessPhiBound.evaluate(jsonoptiEta, betap, doc_ctp, mup, siginvp, sigmaentropyp)
            //print(jsonoptiEta)
            
            //calculated values
            val phi_calculated    = X._1
            val hess_calculated   = X._2._2
            val bound_calculated  = X._3

            //verification values from JSON file
            val bound_test: Double = boundJSON.get(0).toString.toDouble
            val hess_test: DenseMatrix[Double] = DenseMatrix.zeros[Double](74,74) 
            val phi_test: DenseMatrix[Double]  = DenseMatrix.zeros[Double](75,108)
            
            for(i <- 0 until hessianJSON.size()) {
                  var ccc = hessianJSON.get(i).asInstanceOf[JSONArray]     
                  for(j <- 0 until ccc.size()) {
                    hess_test(i,j) = ccc.get(j).toString().toDouble }
            }
                        
            for(i <- 0 until phiJSON.size()) {
                  var ccc = phiJSON.get(i).asInstanceOf[JSONArray]     
                  for(j <- 0 until ccc.size()) {
                    phi_test(i,j) = ccc.get(j).toString().toDouble  }
            }
            
            //compare
            print("\n"+"Bound Calculated		  :"+ X._3) 
            print("\n"+"Bound from JSON		  :"+ bound_test)
                        print("\n"+"Hessian sum of abs diff	: "+sum(abs(hess_calculated - hess_test)))
            print("\n"+"Phi sum of abs diff		: "+sum(abs(phi_calculated - phi_test)))

            /*
            print("\n"+hessianJSON.size()+"x"+hessianJSON.get(0).asInstanceOf[JSONArray].size())
            print("\n"+phiJSON.size()+"x"+phiJSON.get(0).asInstanceOf[JSONArray].size())
            print("\n"+phi_calculated.rows +"x"+ phi_calculated.cols)
            print("\n"+hess_calculated.rows +"x"+ hess_calculated.cols)  
            print("\n"+sum(phi_calculated))
            print("\n"+sum(phi_test))*/

             
          }
      
          def fillInputs() = {     

          //fill inputs from JSON file
          val f  = jsonObject.get("inputs").asInstanceOf[HashMap[String, JSONObject]]
          //println(f.keySet())
           
          var temp = f.get("mu").asInstanceOf[JSONArray] 
           for(i <- 0 until temp.size()) {
            mup(i) = temp.get(i).toString().toDouble
          }
          
          var temp2 = f.get("doc_ct").asInstanceOf[JSONArray] 
           for(i <- 0 until temp2.size()) {
            doc_ctp(i) = temp2.get(i).toString().toDouble
          }
          
          var temp3 = f.get("eta").asInstanceOf[JSONArray] 
           for(i <- 0 until temp3.size()) {
            initialEta(i) = temp3.get(i).toString().toDouble
          }
          
          var temp4 = f.get("beta").asInstanceOf[JSONArray]       
            for(i <- 0 until temp4.size()) {
                var ccc = temp4.get(i).asInstanceOf[JSONArray]     
                for(j <- 0 until ccc.size()) {
                  betap(i,j) = ccc.get(j).toString().toDouble
                }
          }
                
          var temp5 = f.get("siginv").asInstanceOf[JSONArray]       
            for(i <- 0 until temp5.size()) {
                var eee = temp5.get(i).asInstanceOf[JSONArray]     
                for(j <- 0 until eee.size()) {
                  siginvp(i,j) = eee.get(j).toString().toDouble
                }
          }
              
          var temp6 = f.get("sigmaentropy").asInstanceOf[JSONArray]       
             sigmaentropyp = temp6.get(0).toString().toDouble
                     
          }//end fill inputs
  }
  
  //**************************** Spark HDFS [Test] *************************

  def test_sparkhdfs() = {
            import org.apache.spark.{SparkContext, SparkConf}
            import scala.math.random
      
            val conf  = new SparkConf().setAppName("Spark Pi").setMaster("local")
            val spark = new SparkContext(conf)
            val textFile = spark.textFile("hdfs://localhost:9000/whereIsHadoop.txt")
            val counts = textFile.flatMap(line => line.split(" "))
                 .map(word => (word, 1))
                 .reduceByKey(_ + _)
            //counts.saveAsTextFile("hdfs://localhost:9000/CountResults_whereIsHadoop")
            spark.stop()
  }
  
  //**************************** Likelihood [Test] *************************
  def test_likelihood() = {
    println("Maximization of Likelihood Function") 
    
    //JSON import section * start
    try {
      val parser : JSONParser = new JSONParser();
      val reader : FileReader = new FileReader("/Users/aloksingh/Downloads/benchmark2.JSON");
      val jsonObject :JSONObject = (parser.parse(reader)).asInstanceOf[JSONObject]
      
      val outp = jsonObject.get("benchmark").asInstanceOf[HashMap[String, JSONObject]]
      val lhoodJSON : Double = outp.get("lhood").asInstanceOf[JSONArray].get(0).toString.toDouble
       
      /*println("&&&&&&&&&&&&&&&&&&&&& opti eta json normalized &&&&&&&&&&&&&&&&&&&&&")
      println(jsonoptiEtaNorm)
      println("&&&&&&&&&&&&&&&&&&&&& opti eta json normalized &&&&&&&&&&&&&&&&&&&&&")
      */ 
      
      //read inputs from JSON file
      val f  = jsonObject.get("inputs").asInstanceOf[HashMap[String, JSONObject]]
      //println(f.keySet())
       
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
             
      
      val likeliHoodF  = likelihood.lhoodFunction(betap, doc_ctp, mup, siginvp)
      val lbfgs = new LBFGS[DenseVector[Double]](maxIter = 1000000,  m=20) //tolerance = 1E-22, )
      val newEta       = lbfgs.minimize(likeliHoodF, initialEta) 
      val newGradient = likeliHoodF.gradientAt(newEta)
      
      // read gradient @ initial eta from JSON
      var ggg = outp.get("grad").asInstanceOf[JSONArray]
      val jsongradopti : DenseVector[Double] = DenseVector.zeros[Double](ggg.size()) 
      var gradSum = 1.0
      for(i <- 0 until ggg.size()) {
        jsongradopti(i) = (ggg.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
      }
    
      println("Gradient @InitialEta sum of abs diff	: "+sum(abs((jsongradopti - likeliHoodF.gradientAt(initialEta)))))
      println("Likelihood at InitialEta: "+likeliHoodF.valueAt(initialEta))
      println("Likelihood from JSON: "+lhoodJSON)
      
      /*
      // ||||||||||||||||||| 1.1 NORMALIZING newGradient @ *calculated optimal* newEta with this implementation
      var newSum = 1.0
      for(i <- 0 until newGradient.length) {
        newSum = newSum + exp(newGradient(i))
      }
      
      val newgradNorm : DenseVector[Double] = DenseVector.zeros[Double](newGradient.length+1)     
      for(i <- 0 until newGradient.length) {
        newgradNorm(i) = exp(newGradient(i)) / newSum
      }
      newgradNorm(74) = 1.00 / newSum 
      
      // ||||||||||||||||||| 1.2 NORMALIZING JSON final Gradient 
      
      //val outp = jsonObject.get("benchmark").asInstanceOf[HashMap[String, JSONObject]]
      var ggg = outp.get("grad").asInstanceOf[JSONArray]
      val jsongradopti : DenseVector[Double] = DenseVector.zeros[Double](ggg.size()) 
      val jsongradNorm : DenseVector[Double] = DenseVector.zeros[Double](ggg.size()+1)     

      var gradSum = 1.0
      for(i <- 0 until ggg.size()) {
        jsongradopti(i) = (ggg.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
        gradSum = gradSum + exp(jsongradopti(i))
       } 
      
      for(i <- 0 until ggg.size()) {
        jsongradNorm(i) = exp(jsongradopti(i)) / gradSum
      }
      jsongradNorm(74) = 1.0 / gradSum
      
      // ||||||||||||||||||| 2.1 NORMALIZING JSON *Optimal* eta
      var jsonoptiEtaSum = 1.0
      var uuu = outp.get("hpbcpp").asInstanceOf[JSONObject].get("eta").asInstanceOf[JSONObject].get("lambda").asInstanceOf[JSONArray]
      val jsonoptiEta : DenseVector[Double] = DenseVector.zeros[Double](uuu.size()) 
      val jsonoptiEtaNorm : DenseVector[Double] = DenseVector.zeros[Double](uuu.size()+1)     

       for(i <- 0 until uuu.size()) {
        jsonoptiEta(i) = (uuu.get(i).asInstanceOf[JSONArray]).get(0).toString.toDouble
        jsonoptiEtaSum = jsonoptiEtaSum + exp(jsonoptiEta(i))
       }
      
      for(i <- 0 until uuu.size()) {
        jsonoptiEtaNorm(i) = exp(jsonoptiEta(i)) / jsonoptiEtaSum } 
      jsonoptiEtaNorm(74) = 1.0 / jsonoptiEtaSum
      
      // ||||||||||||||||||| 2.2 NORMALIZING new *calculated optimal* eta
      var newEtaSum = 1.0
      val newEtaNorm : DenseVector[Double] = DenseVector.zeros[Double](newEta.length+1)     

       for(i <- 0 until uuu.size()) {
        newEtaSum = newEtaSum + exp(newEta(i))
       }
      
      for(i <- 0 until uuu.size()) {
        newEtaNorm(i) = exp(newEta(i)) / newEtaSum
      }
      newEtaNorm(74) = 1.0 / newEtaSum
      
      
      println("*withOUT* NORMALIZATION")
      println("Optimal Eta sum of abs diff	: "+sum(abs((jsonoptiEta - newEta))))
      println("Final Gradient sum of abs diff	: "+sum(abs((jsongradopti - newGradient))))
      println("Optimal Likelihood from Calculation: "+likeliHoodF.valueAt(newEta))
      println("Optimal Likelihood from JSON: "+lhoodJSON)
       
      println("with NORMALIZATION of L, G, L_json, G_json")
      println("Optimal Eta sum of abs diff	: "+sum(abs((jsonoptiEtaNorm - newEtaNorm))))
      println("Final Gradient sum of abs diff	: "+sum(abs((jsongradNorm - newgradNorm))))
      println("Optimal Likelihood from Calculation: "+likeliHoodF.valueAt(newEta))
      println("Optimal Likelihood from JSON: "+lhoodJSON)
      */ 
      
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