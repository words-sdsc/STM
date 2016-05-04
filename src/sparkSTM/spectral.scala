package sparkSTM


import org.la4j.matrix.sparse.{CRSMatrix => SparseMatrix}
import breeze.linalg.{Axis, DenseMatrix, DenseVector, sum, *, diag} //Matrix, inv, sum, det, cholesky}
import breeze.numerics.sqrt

object spectral {
  
  def gram(mat: SparseMatrix) : DenseMatrix[Double] = {
    val nd = rowSums(mat)
    val indx = nd.findAll { x => x >=2 }
    val matfilter = mat.select(indx.toArray, (0 to mat.columns()-1).toArray)
    val ndfilter = nd(indx).toDenseVector
    val divisor  = ndfilter :* (ndfilter :- 1.0)
    
    //convert mat to densematrix
    val dmat = DenseMatrix.zeros[Double](matfilter.rows(), matfilter.columns())
    for(i <- 0 to matfilter.columns()-1) {
      dmat(::, i) := matfilter.getColumn(i).asInstanceOf[DenseVector[Double]]
    }
    
    /*		val G :DenseMatrix[Double] = (dmat(::, *) :/ divisor) 	//divide each col
      		val F : DenseVector[Double] = sum((dmat(::, *) :/ divisor), Axis._0) */  //Sum down each column (giving a row vector)
    val Htilde = dmat(::, *) :/ sqrt(divisor)        //divide each col
    val Hhat   = diag( sum((dmat(::, *) :/ divisor), Axis._0).toDenseVector)   
    val Q : DenseMatrix[Double] = (Htilde.t * Htilde) - Hhat
  
    Q
    }
 
  def docsToSparseMatrix(documents: List[DenseMatrix[Double]]) : SparseMatrix = {    
    val counts: List[Int]      = documents.map { x => x.cols }  //no of non-zero entries per doc
    val rowPointers: List[Int] = counts.scanLeft(0)(_ + _)      //convert to CRS format
    val colIndices: List[Int]  = documents.map{ x => x(0,::).t.toArray }.flatMap {y => y}.map{t => t.toInt}
    val vals: List[Double]     = documents.map{ x => x(1,::).t.toArray }.flatMap {y => y}
     
     new SparseMatrix(documents.length,colIndices.max,rowPointers.last, vals.toArray, colIndices.toArray, rowPointers.toArray)
     //int rows, int columns, int cardinality, double[] values, int[] columnIndices, int[] rowPointers   
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
  
  //return indices which are zero
  def whichZeros(vec : DenseVector[Double]) : Seq[Int] = {
            var zeroIndices = List[Int]()
            val iter = vec.foreachPair( (i, v)  => {if(v==0) { zeroIndices = i :: zeroIndices } } ) 
            zeroIndices
  }
  
  //return indices that have value >= threshold
  def whichThreshold(vec : DenseVector[Double], thres: Double) : Array[Int] = {
            var whichIndices = List[Int]()
            val iter = vec.foreachPair( (i, v)  => {if(v >= thres) { whichIndices = i :: whichIndices } } ) 
            whichIndices.toArray
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
  
  def fastAnchor(Qbar: DenseMatrix[Double], K: Int, verbose: Boolean): DenseVector[Int] = {
    val basis = DenseVector.zeros[Int](K)
    var rowSquaredSums :DenseVector[Double] = sum( Qbar :* Qbar, Axis._1 ) //col vector
    
    for(i <- 0 to K-1) {
      basis(i) = rowSquaredSums.argmax
      val maxval   = rowSquaredSums.max
      val normalizer = 1/sqrt(maxval)
      Qbar(basis(i), ::) :*= normalizer
      //val rowwisemultiply: DenseMatrix[Double] = Qbar(*,::).map { x => (x :* vec)  } //for each row * with a densevec
      val innerproducts: DenseVector[Double] = Qbar * Qbar(basis(i), ::).t
      val project: DenseMatrix[Double] = innerproducts * Qbar(basis(i), ::)
      basis(0 to i).foreach { x => {project(x,::) := 0.0} } 
      Qbar :-= project
      rowSquaredSums = sum( Qbar :* Qbar, Axis._1 ) //col vector
      basis(0 to i).foreach { x => {rowSquaredSums(x) = 0.0} }
      if(verbose) print(".")
    }
    
    basis
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