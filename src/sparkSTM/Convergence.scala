/**
 * Structural Topic Model by Molly Roberts, Brandon Stewart, Dustin Tingley
 * Website:	http://structuraltopicmodel.com/
 * R Package: https://github.com/bstewart/stm
 * 
 * Scala:
 * @author  Alok Singh (a1singh@ucsd.edu)
 */

package sparkSTM

class Convergence {
  
  var bound : List[Double] = Nil
  var its : Int = 1
  var converged : Boolean = false
  var stopits : Boolean = false
  
}