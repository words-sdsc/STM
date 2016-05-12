package sparkSTM

import org.la4j.matrix.sparse.{CRSMatrix => SparseMatrix}
import breeze.linalg.{Axis, DenseMatrix, DenseVector, sum, *, diag} //Matrix, inv, sum, det, cholesky}
import breeze.numerics.{abs, sqrt, exp}

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
      dmat(::, i) := DenseVector(matfilter.getColumn(i).toDenseVector().toArray()) //.asInstanceOf[DenseVector[Double]]
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
  
  //drop elements from a vector and return remaining vector
  def dropelements(vec: DenseVector[Double], which: Seq[Int]) : DenseVector[Double] = {
           val keep : List[Int] = (0 to vec.length-1 toList) diff which.toList
           var L  = List[Double]()
           for ( i <- keep.reverse ) {
             L = vec(i) :: L 
           }
           DenseVector(L.toArray)      
  }
  
  //fill zeros in a matrix of KxV at columns given by 'whichZeros'
  def refillZeros(K:Int, V:Int, beta0:DenseMatrix[Double], keep: Seq[Int], whichZeros: Seq[Int]) : DenseMatrix[Double] = {
       val betaNew = DenseMatrix.zeros[Double](K, V)
       betaNew(::, keep)       := beta0
       
       //reference: http://lampsvn.epfl.ch/trac/scala/browser/scala/tags/R_2_7_4_final/src/library/scala/Math.scala?view=markup
       betaNew(::, whichZeros) := java.lang.Double.MIN_VALUE 
       
       //divide every col by denominator=row sums
       betaNew(::,*) :/ sum(betaNew(*, ::))
 }
  
  def fastAnchor(Qbar: DenseMatrix[Double], K: Int, verbose: Boolean): List[Int] = {
    var basis = List[Int]()
    var rowSquaredSums :DenseVector[Double] = sum( Qbar :* Qbar, Axis._1 ) //col vector = rowSums
    var indx : Int = 0
    
    for(i <- 0 to K-1) {
      indx                 = rowSquaredSums.argmax
      basis                ::= indx               //adds to 0th position of list, hence reverse later
      
      val maxval           = rowSquaredSums.max
      val normalizer       = 1/sqrt(maxval)
      
      Qbar(indx, ::) :*= normalizer
      
      val innerproducts: DenseVector[Double] = Qbar * Qbar(indx, ::).t
      val project: DenseMatrix[Double]       = innerproducts * Qbar(indx, ::)
      
      project(basis, ::)    := 0.0      
      Qbar                 :-= project
      rowSquaredSums         = sum( Qbar :* Qbar, Axis._1 ) //col vector = rowSums
      rowSquaredSums(basis) := 0.0
      
      if(verbose) print(".")
    }
    
    basis.reverse
  }
  
  def recoverL2(Qbar : DenseMatrix[Double], anchor : List[Int], 
      wprob : DenseVector[Double], verbose : Boolean) : DenseMatrix[Double] = {

    val X         = Qbar(anchor, ::).toDenseMatrix
    val XtX       = X * X.t
    var condprob  = List[DenseVector[Double]]()
    
    def f(i: Int) = {
      
    }
    
    //(0 to Qbar.rows-1).par.foreach { v => f(v) }
    
    for(i <- 0 to Qbar.rows-1) {
      if(anchor.contains(i)) {
        val vec    = DenseVector.zeros[Double](XtX.rows)
        vec(anchor.indexOf(i)) = 1.0
        condprob ::= vec //adds to 0th position of list, hence reverse later
      } else {
        val y      = Qbar(i,::).t
        //println("calling expgrad... for i = " + i)
        condprob ::= spectral.expgrad(X,y,XtX)
      }
      
      if(verbose) { if(i%100==0) println(i + " of " + Qbar.rows) }
    }
    
    if(verbose) println(".exit from recoverL2.")
    
    //take each vector from the list and stack one above another (rbind)
    val weights  = DenseMatrix.vertcat((condprob.reverse).map(_.toDenseMatrix): _*)
    val A        = weights(::,*) :* wprob      //[check] multiply each col of weights by vec wprob
    A.t(::,*) / sum(A, Axis._0).toDenseVector  //sum gives row vector = colSums
    //last line is beta
  }
  
  def expgrad(X : DenseMatrix[Double],y: DenseVector[Double],XtX : DenseMatrix[Double]) : DenseVector[Double] = {
    //note: DenseVector is a col vector; currently y here is passed as a col-vector
    
    val tol     = 1e-7
    val maxIter = 500
    val filler : Double = 1/X.rows
    var alpha   = DenseMatrix.fill(1, X.rows){filler}                            //row-vector
    
    val ytX = new DenseMatrix(1, X.rows , data=(y.t * X.t).t.toArray)            //row-vector
    //[check] is y supposed to be a row vector or a column vector (in R code) ?
    //[check] had to take y.t to make it row vector since ytX is row vector

    var converged = false
    val eta       = 50.0
    var sseOld : Double = java.lang.Double.POSITIVE_INFINITY
    var its       = 1
    
    while(!converged && (its < maxIter)) {
      var grad : DenseMatrix[Double] = ytX - (alpha*XtX)
      val sse  : Double  = sum(grad :* grad)
      grad          = grad :* (eta*2.0)
      val maxderiv  = grad.max
      
      alpha = alpha :* exp(grad - maxderiv)
      //[check] supposed to be element wise multiplication with alpha ?
      
      alpha     = alpha / sum(alpha)
      converged = abs(sqrt(sseOld) - sqrt(sse)) < tol
      sseOld    = sse
      its      += 1
    }
    
    alpha.toDenseVector
    //return(list(par=as.numeric(alpha), its=its, converged=converged,entropy=entropy, log.sse=log(sse)))
  }
  
  def mpinv() = {
    
  }
  
  def tsneAnchor(Q: DenseMatrix[Double]): List[Int] = {
    
    List[Int](0)
  }
   
  def gram_rp() = {
    
  }
  
  
}