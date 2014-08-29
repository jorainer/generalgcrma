normalize.quantiles.exonarray <-
function( x, use.log2=FALSE, verbose=TRUE ){
  if( verbose ) cat("quantile normalising...")
  PMs <- unlist( pmindex( x ), use.names=FALSE )
  PMs <- unique( PMs )
  ## ok, here we go, now i have unique indices of pm probes. can get the expression matrix from those, normalize it and store it back...
  E <- exprs( x )[ PMs, ]
  RN <- rownames(E)
  CN <- colnames(E)
  E <- normalize.quantiles.robust( E, use.log2=use.log2 )
  rownames( E ) <- RN
  colnames( E ) <- CN
  exprs( x )[ PMs, ] <- E
  if( verbose ) cat("finished\n")
  return( x )
}

