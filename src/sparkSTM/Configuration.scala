package sparkSTM
import scala.collection.mutable.{Map}

class Configuration {
  
  var dim   :Map[String, scala.Int] = Map("K"->0, "V"->0, "A"->0, "N"->0)
  var kappa :Map[String, scala.Any] = Map("LDAbeta"->null, "Interactions"->null)
  var init  :Map[String, scala.Any] = Map[String, scala.Any]()
  var convergence  :Map[String, scala.Any] = Map("threshold"->null, "max_iterations"->null)

  var init$verbose = false
  
  init += ("mode"   -> null)
  init += ("nits"   -> null)
  init += ("alpha"  -> null)
  init += ("eta"    -> null)
  init += ("burnin" -> null)
   init += ("s"      -> null)
  init += ("p"      -> null)
  init += ("d-group-size" -> null)
  
}