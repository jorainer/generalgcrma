\name{gcrma.bkg.correct.matrix}
\Rdversion{1.1}
\alias{gcrma.bkg.correct.matrix}
\title{
GCRMA background adjustment
}
\description{
Performs a GCRMA backgorund adjustment on generic input data (i.e. does not require an AffyBatch or CDF environment). Takes a matrix with (optical background) adjusted intensity values and corresponding probe affinities as input and performs a GCRMA NSB-adjustment using functions from the gcrma package.
}
\usage{
gcrma.bkg.correct.matrix(x, affinities, NCprobe, k = 0.5, verbose = TRUE, GSB.adjust = TRUE, cl, mc.cores )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix with (optical background corrected) probe intensity values.
}
  \item{affinities}{
A matrix with the probe affinities. Probe affinities have to be pre-calculated from the probe sequences e.g. using the \code{\link{gcrma.compute.affinities.local.pmonly}} function.
}
  \item{NCprobe}{
Numerical vector with the index of background probes within the input matrix x.
}
  \item{k}{
A tuning factor.
}
  \item{verbose}{
Print progress messages.
}
  \item{GSB.adjust}{
If TRUE, probe effects in specific binding will be adjusted.
}
\item{cl}{ Optional paramter. A cluster object created by the \code{makeCluster} function from the \code{snow} package. If submitted a parallelized version of this function is used. }
\item{mc.cores}{ Optional paramter. Number of cores (CPUs) on which the background adjustment should be performed. Provides, similar to parameter \code{cl}, a significant speed-up of the calculation. Requires that the \code{multicore} package is installed. }
}
\details{
Just a function that can be used to apply GCRMA background adjustment to any (Affymetrix based) microarray.
}
\value{
A matrix with GCRMA background adjusted intensities.
}
\references{

  Rainer J, Lelong J, Bindreither D, Mantinger C, Ploner C, Geley S,
  Kofler R. (2012) Research resource: transcriptional response to
  glucocorticoids in childhood acute lymphoblastic leukemia. Mol Endocrinol. 2012 Jan;26(1):178-93.

  http://bioinfo.i-med.ac.at
}
\author{
Johannes Rainer
}

\seealso{
 \code{\link{gcrma.compute.affinities.local.pmonly}}, \code{\link{gcrma.adjust.nomm}}
}
\examples{
## loading some data
data( "test.sequences" )
data( "affinity.spline.coefs.exon" )
data( "test.data" )

## calculating probe affinities (based on their sequences)
affinity.info <- gcrma.compute.affinities.local.pmonly( sequence=test.sequences, affinity.spline.coefs=affinity.spline.coefs.exon )

## performing a simple optical baclkground adjustment per array (adding +1 to avoid values between 0 and 1)
test.data <- apply( test.data, MARGIN=2, function( x ){ ( x-min( x ) )+1 } )

## just picking 100 random probes in the data and defining them as "non-specific background probes".
## in a real analysis probes in the data should be submitted with the NCprobe that do not measure specific
## signal (e.g. probes not targeting any known transcript)
ncprobes <- sample( 1:nrow( test.data ), 100 )

par( mfrow=c( 1, 2 ) )
plotDensity( log2( test.data[ -ncprobes, ] ), main="specific probes", ylab=expression(log[2]~expression) )
plotDensity( log2( test.data[ ncprobes, ] ), main="background probes", ylab=expression(log[2]~expression) )

## performing the gcrma background adjustment.
affinity.info <- matrix( rep( affinity.info, 2 ), ncol=ncol( test.data ) )
data.bgadjusted <- gcrma.bkg.correct.matrix( test.data, affinities=affinity.info, NCprobe=ncprobes, GSB.adjust=FALSE )
## just adding +1 if any adjusted probe intensity is smaller 1
if( min( data.bgadjusted ) < 1 )
    data.bgadjusted <- data.bgadjusted + 1

plotDensity( log2( cbind( data.bgadjusted[ -ncprobes, ], test.data[ -ncprobes, ] ) ), main="specific probes", xlab=expression(log[2]~expression), col=c( 3, 3, 1, 1 ), lty=c( 1, 1, 2, 2 ) )
legend( "topright", legend=c( "before gcrma", "after gcrma" ), col=c( 1, 3 ), lty=c( 2, 1 ) )


data.bgadjusted.GSB <- gcrma.bkg.correct.matrix( test.data, affinities=affinity.info, NCprobe=ncprobes, GSB.adjust=TRUE )
## just adding +1 if any adjusted probe intensity is smaller 1
if( min( data.bgadjusted.GSB ) < 1 )
    data.bgadjusted.GSB <- data.bgadjusted.GSB + 1

plotDensity( log2( cbind( data.bgadjusted.GSB[ -ncprobes, ], test.data[ -ncprobes, ] ) ), main="specific probes", xlab=expression(log[2]~expression), col=c( 3, 3, 1, 1 ), lty=c( 1, 1, 2, 2 ) )
legend( "topright", legend=c( "before gcrma", "after gcrma" ), col=c( 1, 3 ), lty=c( 2, 1 ) )

}
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
