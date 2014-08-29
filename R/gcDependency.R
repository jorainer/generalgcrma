## some visualisations for the GC dependency...
### added on 7/10/2010 replacement final plotting function of gcBoxplot to an improved one...

if( !isGeneric("gcBoxplot") )
    setGeneric("gcBoxplot", function(x, keep.pattern, exclude.pattern,...)
    standardGeneric("gcBoxplot"))

setMethod( "gcBoxplot", "AffyBatch", function(x, keep.pattern, exclude.pattern,... ){
  #which.probes <- match.arg( which.probes, c( "pm", "mm", "both" ) )
  which.probes <- "pm"
  if( !missing( keep.pattern ) & !missing( exclude.pattern ) ){
    stop( "Either keep.pattern or exclude.pattern should be used, not both!" )
  }
  if( missing( keep.pattern ) & missing( exclude.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes), use.names=FALSE ) )
  }
  if( !missing( keep.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes)[ grep( names( indexProbes( x ) ), pattern=keep.pattern ) ], use.names=FALSE ) )
    if( length( pmIndex )==0 ){
      stop( paste( "No probe sets in the CDF matched the pattern",keep.pattern,"!" ) )
    }
  }
  if( !missing( exclude.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes)[ -grep( names( indexProbes( x ) ), pattern=exclude.pattern ) ], use.names=FALSE ) )
    if( length( pmIndex )==0 ){
      stop( paste( "No probes left after removing all probes from probe sets using the pattern",exclude.pattern,"!" ) )
    }
  }
  library( paste( x@annotation, "probe", sep="" ), character.only=TRUE )
  p <- get( paste( x@annotation, "probe", sep="" ) )
  p <- check.probes( p, x@annotation )

  ## make shure the pm index is unique...
  pmIndex <- unique( pmIndex )
  
  ## have to re-order the sequences, they have to match the ordering of the probes in the CDF file...
  #subIndex <- match( pmIndex, xy2indices( p$x, p$y, cdf=paste( cdfName( x ), "cdf", sep="" ) ) )
  subIndex <- match( pmIndex, xy2indices( p$x, p$y, cdf=paste( x@annotation, "cdf", sep="" ) ) )
  p <- p[ subIndex, ]
  rm( subIndex )

  if( nrow( p )!=length( pmIndex ) ){
    stop("got some internal error...\n")
  }

  ## check if the probe sets in the CDF match with the probe sets in the probe package...
  if( any( pmIndex !=xy2indices( p$x, p$y, cdf=paste( x@annotation, "cdf", sep="" ) ) ) ){
    stop( "probe ordering of the probe package does not match with ordering of the pm probes defined in the CDF file!\n" )
  }
  #doGcBoxplot( intensity( x )[ pmIndex, ], sequences=p$sequence, nt=nt, ylim=ylim )
  ### here my call
  plotNucleoCont( x=exprs(x)[ pmIndex, ], sequences= p$sequence, ...)
  # not sure I need ... there??
} )

##### here the old version
#doGcBoxplot <- function( x, sequences, nt=c( "C", "G" ), ylim=NULL ){
#  AverageInt <- apply( log2( x ), MARGIN=1, mean )
#  Counts <- alphabetFrequency( DNAStringSet( sequences ), baseOnly=TRUE )

#  GCs <- rep( 0, nrow( Counts ) )
#  for( nuc in nt ){
#    GCs <- GCs + Counts[ , nuc ]
#  }
  #GCs <- Counts[ , "C" ] + Counts[ , "G" ]
#  rm( Counts )
#  UniqueCounts <- sort( unique( GCs ) )
#  ## defining the xlab:
#  Xlab <- "Number of "
#  if( length( nt ) > 1 ){
#    NTs <- paste( nt[ 1:( length( nt ) - 1 ) ], collapse=", " )
#    NTs <- paste( c( NTs, nt[ length( nt ) ] ), collapse=" and " )
#    Xlab <- paste( c( Xlab, NTs ), collapse="" )
#  }
#  else{
#    Xlab <- paste( c( Xlab, nt ), collapse="" )
#  }
  
  ## making a boxplot for each of these friends...
#  colfun <- colorRampPalette( brewer.pal( 9, "GnBu" ) )
#  Cols <- colfun( length( UniqueCounts ) )

#  if( is.null( ylim ) ){
#    ylim <- c( min( AverageInt ), max( AverageInt ) )
#  }
#  plot( 3, 3, pch=NA, xlim=c( 1, length( UniqueCounts ) ), ylim=ylim, main="nucleotide count dependent probe intensities", xaxt="n", xlab=Xlab, ylab=expression( log[2]~intensity ) )
#  axis( side=1, at=1:length( UniqueCounts ), label=UniqueCounts )
#  for( i in 1:length( UniqueCounts ) ){
#    boxplot( AverageInt[ GCs==UniqueCounts[ i ] ], add=TRUE, at=i, range=0, col=Cols[ i ] )
#  }
#}

### here the new version
plotNucleoCont <- function (x, sequences, nt = c("C", "G"), log_trans= TRUE, center=FALSE, ylim=NULL, jahcol=FALSE) {
  # x: is an exprs matrix log or natural or a numeric vector of probes intensity
  # sequences: is the corresponding vector of sequence same order as x
  # nt: the nucleotide with possible values {A,T,C,G} (1,2 or 3 as 4 doesn't makes much sense)
  # log_trans: shall we log2 the intensity (default TRUE)
  # jahcol: fancy colramp
  # center: mean probe intensity are centered on the mean(x)
 
  if(log_trans){x<-log2(x)}
  
  # so the easiest way would be to consider x in matrix what ever is the input
  AverageInt<-rowMeans(as.matrix(x))
  # should work even if a numeric vector is inputed
  
  Counts <- alphabetFrequency(DNAStringSet(sequences), baseOnly = TRUE) #fast
  GCs <- rep(0, nrow(Counts))
  ifelse(length(nt)>1, yes=GCs<-GCs + rowSums(Counts[,nt]), no=GCs<-GCs + Counts[,nt])
  rm(Counts)

  if(jahcol){
    colfun<- colorRampPalette(c("red","gold","green"),space="rgb")
  }else{
    colfun <- colorRampPalette(brewer.pal(9, "GnBu"))
  }
  Cols<-colfun(length(unique(GCs)))
  
  # center arround overall mean if center true and set ylab according to it
  if(center){
    M<-mean(x)
    AverageInt<-split(AverageInt-M, GCs )
    ylab=expression(Delta~log[2]~mean~probes~intensity)
  }else{
    AverageInt<-split(AverageInt, GCs )
    expression(log[2]~mean~probes~intensity)
    ylab=expression(log[2]~mean~probes~intensity)
}
  main=paste(paste(nt,collapse="-"),"content bias on probes intensity",sep=" ")
  
  boxplot(AverageInt, main=main, ylim=ylim , xlab =paste(paste(nt,collapse="-"), " count in probes sequences",sep="") , ylab = ylab ,col=Cols,range=0,cex.axis=0.70, varwidth=TRUE, cex.main=0.85) 
}# end of function
######################### done #####################

if( !isGeneric("nucleotidePositionPlot") )
    setGeneric("nucleotidePositionPlot", function(x, keep.pattern, exclude.pattern, ...)
    standardGeneric("nucleotidePositionPlot"))

setMethod( "nucleotidePositionPlot", "AffyBatch", function( x, keep.pattern, exclude.pattern, ... ){
  which.probes <- "pm"
  if( !missing( keep.pattern ) & !missing( exclude.pattern ) ){
    stop( "Either keep.pattern or exclude.pattern should be used, not both!" )
  }
  if( missing( keep.pattern ) & missing( exclude.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes), use.names=FALSE ) )
  }
  if( !missing( keep.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes)[ grep( names( indexProbes( x ) ), pattern=keep.pattern ) ], use.names=FALSE ) )
    if( length( pmIndex )==0 ){
      stop( paste( "No probe sets in the CDF matched the pattern",keep.pattern,"!" ) )
    }
  }
  if( !missing( exclude.pattern ) ){
    pmIndex <- unique( unlist( indexProbes( x, which.probes)[ -grep( names( indexProbes( x ) ), pattern=exclude.pattern ) ], use.names=FALSE ) )
    if( length( pmIndex )==0 ){
      stop( paste( "No probes left after removing all probes from probe sets using the pattern",exclude.pattern,"!" ) )
    }
  }

  library( paste( x@annotation, "probe", sep="" ), character.only=TRUE )
  p <- get( paste( x@annotation, "probe", sep="" ) )
  p <- check.probes( p, x@annotation )

  pmIndex <- unique( pmIndex )
  
  ## have to re-order the sequences, they have to match the ordering of the probes in the CDF file...
  subIndex <- match( pmIndex, xy2indices( p$x, p$y, cdf=paste( x@annotation, "cdf", sep="" ) ) )
  if( length( subIndex )!=length( pmIndex ) ){
    stop( "got some internal error...\n" )
  }
  p <- p[ subIndex, ]
  rm( subIndex )

  if( nrow( p )!=length( pmIndex ) ){
    stop("got some internal error...\n")
  }

  ## check if the probe sets in the CDF match with the probe sets in the probe package...
  if( any( pmIndex !=xy2indices( p$x, p$y, cdf=paste( x@annotation, "cdf", sep="" ) ) ) ){
    stop( "probe ordering of the probe package does not match with ordering of the pm probes defined in the CDF file!\n" )
  }
 # doNucleotidePositionPlot( x=intensity( x )[ pmIndex, ], sequences=p$sequence, FUN=FUN, log.transform=log.transform )
 ## new function call
  plotNucleoPos(x=intensity( x )[ pmIndex, ], sequences=p$sequence , ...)
} )
######## old version
#doNucleotidePositionPlot <- function( x, sequences, FUN=mean, log.transform=TRUE ){
#  if( log.transform ){
#    x <- log2( x )
#  }
#  AverageInt <- apply( x, MARGIN=1, FUN=FUN )
  ## making a matrix with nucleotides... each row is one probe.
#  SeqMat <- t( sapply( sequences, function( z ){ unlist( strsplit( z, split="" ), use.names=FALSE ) }, USE.NAMES=FALSE ) )
#  NucAverage <- matrix( nrow=4, ncol=ncol( SeqMat ) )
#  rownames( NucAverage ) <- c( "A", "C", "G", "T" )
#  for( nucleotide in rownames( NucAverage ) ){
#    for( i in 1:ncol( NucAverage ) ){
      #cat( nucleotide,":",i, "of 25\n" )
#      NucAverage[ nucleotide, i ] <- mean( AverageInt[ SeqMat[ , i ]==nucleotide ], na.rm=TRUE )
#    }
#  }
  ## finally generating the plot...
#  Cols <- brewer.pal( 4, "Set1" )
#  names( Cols ) <- c( "A", "C", "G", "T" )
#  plot( 3, 3, pch=NA, xlim=c( 1, ncol( NucAverage ) ), ylim=c( min( NucAverage, na.rm=TRUE ), max( NucAverage, na.rm=TRUE ) ), xlab="nucleotide position", ylab="average intensity", main="Nucleotide position dependent average intensity" )
#  abline( h=0, col="grey", lty=2 )
#  for( nuc in rownames( NucAverage ) ){
#    points( x=1:ncol( NucAverage ), y=NucAverage[ nuc, ], pch=nuc, col=Cols[ nuc ], type="b" )
#  }
#}

### here new version (no for loop faste averageInt, allow plot on numeric vector and centered)

plotNucleoPos<- function (x,
                          sequences,
                          log_trans = TRUE,
                          center=FALSE) 
{
  if (log_trans) {x <- log2(x)}

  AverageInt<-rowMeans(as.matrix(x))
        
  ## get a matrix row nb of seq , col letters@pos from splited sequence (25) 
  SeqMat <- t(sapply(sequences, function(z) {unlist(strsplit(z, split = ""), use.names = FALSE)}, USE.NAMES = FALSE))# this is brillant thanks to jo
  ##set the mat to receive mean probe int letter@pos
  NucAverage <- matrix(nrow = 4, ncol = ncol(SeqMat))
  rownames(NucAverage) <- c("A", "T", "G", "C")
  M<-mean(x) # we need that guy in both case (to center or to diplay the Mean if no center)
  if(center){
    NucAverage<-t(sapply(rownames(NucAverage),function(base) {
      sapply(1:ncol(NucAverage),function(pos) mean(AverageInt[SeqMat[,pos]==base])-M)}))
    ylab<-expression(Delta~log[2]~mean~probes~intensity)
  }else{
    NucAverage<-t(sapply(rownames(NucAverage),function(base) {
      sapply(1:ncol(NucAverage),function(pos) mean(AverageInt[SeqMat[,pos]==base]))}))
    ylab<-expression(log[2]~mean~probes~intensity)
  }
    
  ##plots
  Cols <- brewer.pal(4, "Set1")
  bases<-c("A", "T", "G", "C")
  plot(NucAverage["A",], type="n", xlim = c(1, ncol(NucAverage)), ylim = c(min(NucAverage),max(NucAverage)), xlab = "nucleotide position in probe sequence", ylab=ylab, main = "Nucleotide position bias", cex.axis=0.7,cex.main=0.85)
  abline(h = ifelse(center,yes=0,no=M), col = "grey", lty = 2)
  ## use tmp because sapply output some empty list... not always... strange
  tmp<-sapply(1:length(bases),function(x) {points(NucAverage[x, ],pch=bases[x],col=Cols[x],type="b" )})
}# end of function

