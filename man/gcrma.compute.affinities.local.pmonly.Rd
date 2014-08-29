\name{gcrma.compute.affinities.local.pmonly}
\Rdversion{1.1}
\alias{gcrma.compute.affinities.local.pmonly}
\title{
Probe Affinity computation
}
\description{
Compute probe affinities based on their sequence.
}
\usage{
gcrma.compute.affinities.local.pmonly(sequence, verbose = TRUE, affinity.spline.coefs = NULL, pmonly = TRUE, mc.cores )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{
  A character vector with the sequences of the individual probes.
}
  \item{verbose}{
  Plot messages about the progress of the function.
}
  \item{affinity.spline.coefs}{
  Affinity spline coefficients (pre-calculated).
}
  \item{pmonly}{
  Currently only TRUE is supported.
}
\item{mc.cores}{
	Number of cores (CPUs) on which the affinity computations should be performed. Requires that the multicore package is installed.
}
}
\details{
See documentation of the \code{compute.affinities.local} function in the gcrma package.
}
\value{
	A matrix with the same number of columns than affinity.spline.coefs (or only one column if affinity.spline.coefs is a vector) containing the affinity scores for the individual probes.
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

\examples{
data( "test.sequences" )
data( "affinity.spline.coefs.exon" )
test.sequences[ 1:4 ]

affinity.info <- gcrma.compute.affinities.local.pmonly( sequence=test.sequences, affinity.spline.coefs=affinity.spline.coefs.exon )

affinity.info[ 1:4 ]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
