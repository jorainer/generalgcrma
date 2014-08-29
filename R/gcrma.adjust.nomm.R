#####
# 08.10.2012: returns also aggregated residual information per probe set.
# 20.07.2012: added option normalize.method that could allow us to run a different normalization method.
# 24.08.2011: added option mc.cores, allows parallel computing using the multicore package.
# 17.02.2011: added option summarize.method to allow either rma or plm fit.
# 25.06.2009
# Performs the gcrma background adjustment for Affy chips without mm probes.
gcrma.adjust.nomm <- function( x,
                              optical.correct=TRUE,
                              verbose=TRUE,
                              affinity.spline.coefs,
                              NCprobe,
                              k=0.5,
                              GSB.adjust=TRUE,
                              NC.pattern="bg",
                              normalize.quantiles=TRUE,
                              normalize.method="quantiles",
                              normalize.log2=FALSE,
                              summarize=FALSE,
                              summarize.method="rma",
                              cl,
                              mc.cores,
                              ret.residual.aggr=FALSE,
                              residual.aggr.fun=function( z ){ mean( abs( z ) ) }
                              ){
  LocalAffinities <- FALSE
  summarize.method <- match.arg( summarize.method, c( "rma", "plm" ) )
  if( ret.residual.aggr ){
      ## have to use PLM, since RMA does not return residuals
      summarize.method="plm"
  }
  ## check for possible values of normalize.method:
  normalize.method <- match.arg( normalize.method, c( "quantiles", "subset.quantiles" ) )

  if( missing( affinity.spline.coefs ) ){
    LocalAffinities <- TRUE
    #cat("affinity.spline.coefs have to be submitted!")
  }
  if( optical.correct ){
#    x <- bg.adjust.optical( x , verbose=verbose )
    ## optical bg adjustment is fast anyways... will not run this in parallel mode...
    x <- gcrma.bg.adjust.optical( x, verbose=verbose )
  }
  if( verbose ) cat( "verifying that probe ordering matches probe sequence ordering\n" )

  # extract the perfect match probe intensities (NOTE: must also contain NSB probes!)
  pmIndex <- unlist( indexProbes( x, "pm") )

  ## need the sequences...
  library( paste( x@annotation, "probe", sep="" ), character.only=TRUE )
  data( list=paste( x@annotation, "probe", sep="" ) ) # since R-3.0 we have to load the data first.
  p <- get( paste( x@annotation, "probe", sep="" ) )
  #rm( list=paste( paste( x@annotation, "probe", sep="" ) ) ) # just to clean a little memory
  p <- check.probes( p, cdfName( x ) )

  ## have to re-order the sequences, they have to match the ordering of the probes in the CDF file...
  subIndex <- match( pmIndex, xy2indices( p$x, p$y, cdf=paste( cdfName( x ), "cdf", sep="" ) ) )
  p <- p[ subIndex, ]
  rm( subIndex )

  ## check if the probe sets in the CDF match with the probe sets in the probe package...
  if( any( pmIndex !=xy2indices( p$x, p$y, cdf=paste( cdfName( x ), "cdf", sep="" ) ) ) ){
    stop( "probe ordering of the probe package does not match with ordering of the pm probes defined in the CDF file!\n" )
  }

  ## Negative control probes...
  if( missing( NCprobe ) ){
    ## try to determine NSB probes from the names of the probe sets...
    #NCprobe <- unique( unlist( pmindex( x )[ grep( names( pmindex( x ) ), pattern=NC.pattern ) ] ) ) ## are the indices for the intensity( x )
    if( verbose )
      cat( paste("Defining NCprobes using a pattern search for", NC.pattern, "in probe set names.\n") )
    NCprobe <- grep( names( pmIndex ), pattern=NC.pattern )
    if( length( NCprobe ) < 1 ){
      stop( paste( "No NCprobes could be found!\n" ) )
    }
  }

  ## get the intensities for the pm probes...
  PM <- intensity( x )[ pmIndex, ]

  ## calculating the affinity spline coefs from within the data...
  if( LocalAffinities ){
    if( verbose ) cat( "estimating affinity spline coefficients from the expression intensities and probe sequences of the background probes..." )
    affinity.spline.coefs <- base.profiles( PM[ NCprobe, ], p$sequence[ NCprobe ] )
    if( verbose) cat( "finished.\n" )
  }

  ## finally calculating the probe affinities... can be shure that the ordering of the probes sequences matches the
  ## ordering of the probes in the PM matrix.
  affinity.info <- gcrma.compute.affinities.local.pmonly( sequence=p$sequence, verbose=verbose, affinity.spline.coefs=affinity.spline.coefs, mc.cores=mc.cores )

  ## just checking that the number of columns of the affinity.info matches the number of columns in the AffyBatch...
  if( ncol( affinity.info )==1 ){
    ## would be faster with affinity.info <- matrix( rep( affinity.info, length( sampleNames( x ) ) ), ncol=lenght(sampleNames( x )) )
    tmp <- affinity.info
    affinity.info <- matrix( ncol=length( sampleNames( x ) ), nrow=nrow( tmp ) )
    for( i in 1:ncol( affinity.info ) ){
      affinity.info[ , i ] <- tmp[ , 1 ]
    }
    rm( tmp )
  }
  else{
    if( ncol( affinity.info )!=length( sampleNames( x ) ) ){
      stop( "nr of columns of the AffyBatch and the affinity.spline.coefs do not match!\n" )
    }
  }

#  rm( pmIndex )
  rm( p )

  # performing the NSB background adjustment
  # NOTE: NCprobe have to match within the data matrix!
  PM <- gcrma.bkg.correct.matrix( x=PM, affinities=affinity.info, NCprobe=NCprobe, k=k, verbose=verbose, GSB.adjust=GSB.adjust, cl=cl, mc.cores=mc.cores )

  ## just adding +1 to all intensities...
  if( min( PM ) <= 1 ){
    if( verbose ) cat("adding +1 to all intensities to avoid negative intensities after log2 transformation\n")
    PM <- PM + 1
  }
  pm( x ) <- PM
  rm( PM )

  if( normalize.quantiles ){
      if( normalize.method=="subset.quantiles" ){
          ## do subset quantiles normalisation as implemented in the SQN package.
          warning( "Subset quantiles normalization not implemented yet! Using quantiles normalization instead." )
          normalize.method <- "quantiles"
      }
      if( normalize.method=="quantiles" ){
          ## quantile normalising...
          x <- normalize.quantiles.exonarray( x, use.log2=normalize.log2, verbose=verbose )
      }
  }

  if( summarize ){
    if( summarize.method=="rma" ){
      if( verbose ){
        cat( "summarizing probe intensities using rma (median polish).\n" )
      }
      x <- rma( x, normalize=FALSE, background=FALSE, verbose=verbose )
    }
    if( summarize.method=="plm" ){
      if( verbose ){
        cat( "summarizing probe intensities using fitPLM.\n" )
      }
      require( affyPLM )
      #subset <- NULL
      #psnames <- unique( names( pmIndex ) )
      #subset <- psnames[ -grep( psnames, pattern=NC.pattern ) ]
      #PLM <- fitPLM( x, output.param=list( residuals=ret.residual.aggr, weights=FALSE ), background=FALSE, normalize=FALSE, subset=subset )
      PLM <- fitPLM( x, output.param=list( residuals=ret.residual.aggr, weights=FALSE ), background=FALSE, normalize=FALSE )
      cat( "finished fitting with fitPLM\n" )
      if( ret.residual.aggr ){
          ## that's now a little more complicated:
          ## 1) extract residuals and run the aggregate function to aggregate residuals per probe set
          ## 2) build an extended version of ExpressionSet and return that.
          cat( "aggregating residuals per probe set..." )
          x <- newResExpressionSet( PLMset2exprSet( PLM ) )
          Resids <- residuals( PLM )$PM.resid
          Resids <- Resids[ rownames( Resids ) %in% featureNames( x ), ]
          Resids <- as.matrix( aggregate( Resids, by=list( probeset=factor( rownames( Resids ) ) ), FUN=residual.aggr.fun )[ , -1 ] )
          rownames( Resids ) <- featureNames( x )
          x@residuals <- Resids
          cat( "done\n" )
      }else{
          ## transform the PLMset into a ExpressionSet
          cat( "transforming the PLMset..." )
          cdfname <- x@cdfName
          x <- PLMset2exprSet( PLM )
          x@annotation <- cdfname
      }
    }
  }
  return( x )
}

