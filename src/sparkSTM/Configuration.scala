package sparkSTM
import scala.collection.mutable.{Map}

class Configuration {
  
  var dim   :Map[String, scala.Any] = Map("K"->null, "V"->null, "A"->null, "N"->null)
  var kappa :Map[String, scala.Any] = Map("LDAbeta"->null, "Interactions"->null)
  var init  :Map[String, scala.Any] = Map[String, scala.Any]()
  
  init += ("mode"   -> null)
  init += ("nits"   -> null)
  init += ("alpha"  -> null)
  init += ("eta"    -> null)
  init += ("burnin" -> null)
  init += ("verbose"-> null)
  init += ("s"      -> null)
  init += ("p"      -> null)
  init += ("d-group-size" -> null)
  
}