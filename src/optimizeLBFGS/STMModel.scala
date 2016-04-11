package optimizeLBFGS

import org.apache.spark.rdd.RDD
import org.apache.spark.mllib.linalg.{Vector}
import breeze.linalg.{DenseMatrix, DenseVector}

class STMModel {
  //GG = set of global parameters
  
  
  //document storage
  private var docs: RDD[(Long, Vector)]= null
  
  /* ********************* run STM over a set of documents ********************* */
  def runSTM(setOfDocs: RDD[(Long, Vector)]): STMModel = {
    
    //get copy of global parameters
    
    //broadcast large data structures to each node if they are needed by every doc
    
    //partition set of documents
    val metricsFromPartitions: RDD[(DenseMatrix[Double], List[DenseVector[Double]])] =   
      setOfDocs.mapPartitions 
      {
        docs =>
                val size             = 5                       //random for now
                // TP = set of accumulators for this partition
                val partition_accum1 = DenseMatrix.zeros[Double](size,size)
                var partition_accum2 = List[DenseVector[Double]]()
                
                docs.foreach 
                {
                    case(id: Long, termCounts: Vector) =>
                      /* send this document for STM processing,
                         DD = set of results from STM for this single document */
                      val i = 1
                    
                    
                    // update partition accumulators TP with DD
                    
                }

                // return iterator (each tuple is a set of accumulators from one partition)
                Iterator((partition_accum1, partition_accum2))
        
      }//end mapPartitions
    
    /*each tuple in 'metricsFromPartitions' is from one partition
      PP = aggregation of results from each partition in 'metricsFromPartitions'
      */
      
    
    
    //update global parameters GG using these aggregates PP
    
    
    //unpersist broadcasted variables
    
    
    this
    
  } //end runSTM
  
  
  /* *********************initialize global parameters ********************* */
  def initialize() = {
    
  }//end initialize
  
}