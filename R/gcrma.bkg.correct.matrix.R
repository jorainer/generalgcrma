#######
# 24.08.2011: option mc.cores: perform parallel computation using the multicore package.
# 25.06.2009
# the function takes a matrix with expression intensities (optical background corrected) and affinities for
# the coresponding probe sequences as input and performs a gcrma background correction using functions from the
# gcrma package (bg.adjust.affinities).
# the matrix MUST contain also values for non-specific probes and the indices of these probes within the input
# matrix x and the affinities matrix has to be submitted.
# ordering of rows in x and affinities has to be the same!
# Note: affinities have to be pre-calculated, rows in affinities have to match rows in the intensity matrix x!
gcrma.bkg.correct.matrix <-
function( x,
         affinities,
         NCprobe,
         k=0.5,
         verbose=TRUE,
         GSB.adjust=TRUE,
         cl,
         mc.cores
         ){
  ## just check that x and affinities have the same number of rows...
  if( missing( NCprobe ) ){
    stop( "indices of non-specific probes have to be submitted using the NCprobe parameter!\n" )
  }
  if( nrow( x )!=nrow( affinities ) ){
    stop( "x and affinities must have the same number of rows" )
  }
  if( ncol( x )!=ncol( affinities ) ){
    stop( "x and affinities must have the same number of columns" )
  }
  index.affinities=1:nrow( x )     # actuall i don't know for what this index is good for...
  ## checking if NCprobe indices are within nrow(x)
  bm <- ! NCprobe %in% index.affinities
  if( any( bm ) ){
    #NCprobe <- NCprobe[ !bm ]
    stop( paste( sum(bm), "background probe indices are larger than the actual data matrix!\n" ) )
  }
  fit1 <- NULL
  if( GSB.adjust ){
    set.seed(1)
    sample.size <- 50000
    if( sample.size > nrow( x ) ){
      sample.size <- floor( nrow( x ) / 10 )
    }
    Subset <- sample( 1:length( as.matrix( x )[ index.affinities, ] ), sample.size)
    Y <- log2( x )[ index.affinities , ][ Subset ]
    #Subset <- (Subset - 1)%%nrow(as.matrix( x )[index.affinities, ]) + 1   # i think i would only need this if affinities was a vector...
    X <- affinities[ index.affinities , ][ Subset ]
    fit1 <- lm(Y~X)$coef
  }
  ## just to be on the save side... keep rownames and colnames...
  CN <- colnames( x )
  #RN <- rownames( x )
  if( verbose ) cat("Adjusting for non-specific binding")
  if( !missing( cl ) | !missing( mc.cores ) ){
    if( verbose ) cat(" in parallel processing mode ")
    ## building a list with values and affinities...
    List <- vector( "list", ncol( x ) )
    for( i in 1:ncol( x ) ){
      List[[i]] <- list( values=x[ , i ], affinities=affinities[ , i ] )
    }
    if( !missing( cl ) ){
      ## use snow
      x <- parSapply( cl, List, .performBgAdjustAffinities, NCprobe=NCprobe, k=k, fast=FALSE, nomm=TRUE, verbose=verbose, GSB.adjust=GSB.adjust, fit1=fit1 )
    colnames( x ) <- CN
    }
    if( !missing( mc.cores ) ){
      ## use multicore/parallel
      tmp <- mclapply( List, FUN=.performBgAdjustAffinities, NCprobe=NCprobe, k=k, fast=FALSE, nomm=TRUE, verbose=verbose, GSB.adjust=GSB.adjust, fit1=fit1, mc.cores=mc.cores, mc.set.seed=FALSE )
      x <- do.call( cbind, tmp )
      rm( tmp )
      colnames( x ) <- CN
    }
  }
  else{
    for( i in 1:ncol( x ) ){
      if( verbose ) cat(".")
      ## ok, that's where we call the function from the gcrma package.
      x[,i] <- bg.adjust.affinities( pms=x[ , i ],
                                    ncs=x[ NCprobe, i ],
                                    apm=affinities[ , i ],
                                    anc=affinities[ NCprobe , i ],
                                    index.affinities=index.affinities,
                                    k=k,
                                    fast=FALSE,
                                    nomm=TRUE
                                    )
      if( GSB.adjust ){
        ##      if( verbose ) cat("Adjusting GSB (gene specific binding)")
        x[,i] <- GSB.adj( Yin=x[,i], subset=index.affinities, aff=affinities[, i],fit1=fit1,k=k)
      }
    }
  }


  if( verbose ) cat( "finished\n" )
  return( x )
}

.performBgAdjustAffinities <- function( x, NCprobe, k, fast, nomm, verbose, GSB.adjust, fit1 ){
  require( "gcrma", quietly=TRUE )
  if( verbose ) cat(".")
  z <- bg.adjust.affinities( pms=x[["values"]],
                            ncs=x[["values"]][ NCprobe ],
                            apm=x[["affinities"]],
                            anc=x[["affinities"]][ NCprobe ],
                            k=k,
                            fast=fast,
                            nomm=nomm,
                            index.affinities=1:length( x[["values"]] )
                            )
  if( GSB.adjust ){
    z <- GSB.adj( Yin=x[["values"]],
                    subset=1:length( x[["values"]] ),
                    aff=x[["affinities"]],
                    fit1=fit1,
                    k=k
                    )
  }
  gc()
    return( z )
}


