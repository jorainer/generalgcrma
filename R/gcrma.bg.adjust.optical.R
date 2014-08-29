## optical bg adjustment... similar to bg.adjust.optical, with the major difference, that all probes on the array are used...
gcrma.bg.adjust.optical <- function( abatch, minimum=1, which.probes="all", verbose=TRUE, cl ){
  which.probes <- match.arg( which.probes, c( "all", "both", "pm", "mm" ) )
  if( which.probes=="all" ){
    Index <- 1:nrow( intensity( abatch ) )
  }
  else{
    Index <- unlist( indexProbes( abatch, which.probes ), use.names=FALSE )
  }
  Index <- Index[ !is.na( Index ) ]
  exprs( abatch )[ Index, ] <- gcrma.bg.adjust.optical.matrix( exprs( abatch )[ Index, ], minimum=minimum, verbose=verbose, cl=cl )
  return( abatch )
} 


gcrma.bg.adjust.optical.matrix <- function( x, minimum=1, verbose=TRUE, cl ){
  if (verbose) cat("Adjusting for optical effect")
  CN <- colnames( x )
  if( missing( cl ) ){
    for( i in 1:ncol( x ) ){
      if( verbose ) cat(".")
      x[ , i ] <- x[ , i ] - min( x[ , i ], na.rm=TRUE ) + minimum
    }
  }
  else{
    ## parallel mode!
    if( verbose ) cat(" in parallel mode...")
    x <- parApply( cl=cl, X=x, MARGIN=2, FUN=function( z ){
      results <- z - min( z, na.rm=TRUE ) + minimum
      return( results )
    } )
    return( x )
  }
  if( verbose ) cat( "finished\n" )
#  colnames( x ) <- CN
  return( x )
}







