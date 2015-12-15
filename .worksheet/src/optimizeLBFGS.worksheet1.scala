package optimizeLBFGS

import breeze.linalg._
import breeze.optimize._

object worksheet1 {;import org.scalaide.worksheet.runtime.library.WorksheetSupport._; def main(args: Array[String])=$execute{;$skip(105); 
		val a = 5.0;System.out.println("""a  : Double = """ + $show(a ));$skip(52); 
		val n : DenseVector[Double]= DenseVector(1.2,3.2);System.out.println("""n  : breeze.linalg.DenseVector[Double] = """ + $show(n ))}
		
/*	val a = DenseMatrix.rand(2,2)
	val b = DenseVector.rand(5)
	
		val dm = DenseMatrix((1.0,2.0,3.0), (4.0,5.0,6.0))
    val re = dm(::,*) :* DenseVector(10.0, 100.0)
    
    val sumRows = sum(dm(::,*))
		sumRows.rows
		sumRows.cols
		
		val sumVector = sumRows.toDenseVector
    sumVector.length
		val other = sumRows.t.toDenseVector
		other.length
 
		val eg = DenseVector(4.3,4.5)
		val anotherV = DenseVector(1.0, 6.5, 4.6)
		val k = anotherV.asDenseMatrix
		k.rows
		k.cols
		
    val a = eg.toArray.toList
    val b = List(9.0)
    val c = a ::: b
    c.toArray*/
}
