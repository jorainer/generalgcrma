\name{gcrma.bg.adjust.optical}
\Rdversion{1.1}
\alias{gcrma.bg.adjust.optical}
\alias{gcrma.bg.adjust.optical.matrix}
\title{
Optical background correction
}
\description{
Adjusts for optical background by subtracting the minimum intensity per microarray (column).
}
\usage{
gcrma.bg.adjust.optical( abatch, minimum=1, which.probes="all", verbose = TRUE, cl)
gcrma.bg.adjust.optical.matrix( x, minimum=1, verbose = TRUE, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{abatch}{ AffyBatch object.}
\item{x}{ A numerical matrix with intensities. }
\item{minimum}{ Minimal intensity value that is added to each bg adjusted intensity (defaults to +1 and thus avoids adjusted values between 0 and 1).}
\item{which.probes}{ Which probes should be bg adjusted. Allowed values are \code{all} (use all probes/intensities on the array), \code{both} (all probes on the array that are either pm or mm probes), \code{pm} (only the pm probes), \code{mm} (only the mm probes).}
  \item{verbose}{
Print progress messages.
}
\item{cl}{ Optional paramter. A cluster object created by the \code{makeCluster} function from the \code{snow} package. If submitted a parallelized version of this function is used. }
}
\details{

}
\value{
A AffyBatch object with (optical) background adjusted intensities.
A matrix with (optical) background adjusted intensities.
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
data( "test.data" )

## performing a simple optical baclkground adjustment per array (adding +1 to avoid values between 0 and 1)
test.data <- gcrma.bg.adjust.optical.matrix( test.data )


}
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
