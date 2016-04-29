package sparkSTM

import org.la4j.matrix.sparse.{CCSMatrix => SparseMatrix}
import breeze.linalg.{DenseMatrix, DenseVector, sum, *} //Matrix, diag, inv, sum, det, cholesky}

object spectral {
  
  def gram(mat: SparseMatrix) : DenseMatrix[Double] = {
    
    DenseMatrix.rand[Double](2, 2) //dummy
  }
  
  def docsToSparseMatrix(documents: List[DenseMatrix[Double]]) : SparseMatrix = {    
     val docsijv : Tuple3[Array[Int], Array[Int], Array[Double]] = null // get ijv from documents
     // docsijv = (Array[Int] colPtrs, Array[Int] rowIndices, Array[Int] values)
     val numRows = 100
     val numCols = 100
     
     //rows = documents, cols = word indices, values = counts
     //int rows, int columns, int cardinality, double[] values, int[] rowIndices, int[] columnPointers
     new SparseMatrix(numRows, numCols)   
  }
  
  def colSums(mat : SparseMatrix) : DenseVector[Double] = {
            //val iter = mat.iteratorOrNonZeroColumns()
            val colSum = DenseVector.zeros[Double](mat.columns())
            var j =0
            while (j < mat.columns()) {
               //val i = iter.next()
               colSum(j) = mat.getColumn(j).sum()
               j=j+1
            }
     colSum
  }
  
  def rowSums(mat : SparseMatrix) : DenseVector[Double] = {
            val rowSum = DenseVector.zeros[Double](mat.rows())
            var i = 0
            while (i < mat.rows()) { 
              rowSum(i) = mat.getRow(i).sum()
              i=i+1
            }
      rowSum
  }
  
  def whichZeros(vec : DenseVector[Double]) : Seq[Int] = {
            var zeroIndices = List[Int]()
            val iter = vec.foreachPair( (i, v)  => {if(v==0) { zeroIndices = i :: zeroIndices } } ) 
            zeroIndices
  }
  
  def dropelements(vec: DenseVector[Double], which: Seq[Int]) : DenseVector[Double] = {
           val keep : List[Int] = (0 to vec.length-1 toList) diff which.toList
           var L  = List[Double]()
           for ( i <- keep.reverse ) {
             L = vec(i) :: L 
           }
           DenseVector(L.toArray)      
  }
  
  def refillZeros(K:Int, V:Int, beta0:DenseMatrix[Double], keep: Seq[Int], whichZeros: Seq[Int]) : DenseMatrix[Double] = {
       val betaNew = DenseMatrix.zeros[Double](K, V)
       betaNew(::, keep)       := beta0
       betaNew(::, whichZeros) := 5.960465E-8 //reference: https://issues.scala-lang.org/browse/SI-3791
       //divide every col by denominator=row sums
       betaNew(::,*) :/ sum(betaNew(*, ::))
 }
  
  def gram_rp() = {
    
  }
  
  def fastAnchor(Q: DenseMatrix[Double], K: Int, verbose: Boolean): DenseVector[Int] = {
    
    DenseVector[Int](0)
  }
  
  def recoverL2(Q : DenseMatrix[Double], anchor : DenseVector[Int], 
      wprob : DenseVector[Double], verbose : Boolean) : DenseMatrix[Double] = {
    
    DenseMatrix.zeros[Double](2,2)
  }
  
  def expgrad() = {
    
  }
  
  def mpinv() = {
    
  }
  
  def tsneAnchor(Q: DenseMatrix[Double]): DenseVector[Int] = {
    
    DenseVector[Int](0)
  }
  
  
}