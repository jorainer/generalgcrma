########
# 24.08.2011: option mc.cores: perform parallel processing...
# 25.06.2009.
# the following is a modification of "compute.affinities.local" from the gcrma package (file base.profiles.2.R),
# designed to bypass CDF and ProbePackage configuration. The function works for PM only chips (e.g. Exon or Tiling
# microarrays).
# The function takes a vector with (PM) sequences and affinity.spline.coefs as input and returns a matrix with
# the affinities for the probes (the matrix has the same number of columns than affinity.spline.coefs).
# the function will return a list with two matrices for pmonly=FALSE...
gcrma.compute.affinities.local.pmonly <-
function( sequence,
         verbose=TRUE,
         affinity.spline.coefs=NULL,
         pmonly=TRUE,
         mc.cores
         ){
  #require( splines, quiet=TRUE )
  if( !pmonly ){
    stop("pmonly=FALSE is not yet implemented...")
  }
  if(is.null(affinity.spline.coefs)){
    stop( "affinity.spline.coefs have to be submitted! please calculate them with the compute.affinity.coef function or submit pre-caluclated affinities.\n" )
  }
  #require( splines, quietly=TRUE )
  if( is.null( ncol( affinity.spline.coefs ) ) ){
    nCOL <- 1
  }
  else{
    nCOL <- ncol( affinity.spline.coefs )
  }
  affinity.basis.matrix <- ns(1:25,df=nrow(as.matrix(affinity.spline.coefs))/3)
  
  A13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[1:5])
  T13 <- 0
  C13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[6:10])
  G13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[11:15])

  if(verbose) cat("Computing affinities")
  
  APM=matrix(NA,length( sequence ),nCOL )
  if( !pmonly ){
    AMM=matrix(NA,length( sequence ),nCOL )
  }
  for( K in 1:ncol(APM)){
    if(verbose) cat(".")
    apm <- vector("numeric",length( sequence ))
    if( !pmonly ){
      amm <- vector("numeric",length( sequence ))
    }
    if( !missing( mc.cores ) ){
      if( nCOL > 1 ){
        aff.spline.coefs <- affinity.spline.coefs[,K]
      }
      else{
        aff.spline.coefs <- affinity.spline.coefs
      }
      apm <- unlist( mclapply( as.list( sequence ), FUN=.doCalculation, mc.cores=mc.cores, mc.set.seed=FALSE, aff.basis.matrix=affinity.basis.matrix, aff.spline.coefs=aff.spline.coefs ) )
    }
    else{
      for(i in seq(along=apm)) {
        charMtrx <- .Call("gcrma_getSeq", sequence[i],
                          PACKAGE="gcrma")
        A <- cbind(charMtrx[1,] %*% affinity.basis.matrix,
                   charMtrx[2,] %*% affinity.basis.matrix,
                   charMtrx[3,] %*% affinity.basis.matrix)
        
        if( nCOL > 1 ){
          apm[i] <- A %*% affinity.spline.coefs[,K]
        }
        else{
          apm[i] <- A %*% affinity.spline.coefs
        }
#      if (charMtrx[1,13] == 1) {
#        amm[i] <- apm[i] + T13 - A13
#      }
#      else {
#        if (charMtrx[4,13] == 1) {
#          amm[i] <- apm[i] + A13 - T13
#        }
#        else{
#          if (charMtrx[3,13]) {
#            amm[i] <- apm[i] + C13 - G13
#          }
#          else {
#            amm[i] <- apm[i] + G13 - C13
#          }
#        }
#      }
      }
    }
    APM[,K]=apm
    #;AMM[,K]=amm
  }
  if(verbose) cat("Done.\n")
  rownames( APM ) <- names( sequence )
  return( APM )
}

.doCalculation <- function( x, aff.basis.matrix, aff.spline.coefs ){
  charMtrx <- .Call("gcrma_getSeq", x,
                    PACKAGE="gcrma")
  A <- cbind(charMtrx[1,] %*% aff.basis.matrix,
             charMtrx[2,] %*% aff.basis.matrix,
             charMtrx[3,] %*% aff.basis.matrix)
  return( A %*% aff.spline.coefs )      
}

