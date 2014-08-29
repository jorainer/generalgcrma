\name{affinity.spline.coefs.exon}
\Rdversion{1.1}
\alias{affinity.spline.coefs.exon}
\docType{data}
\title{
Affinity spline coefficients for Exon microarrays
}
\description{
Pre-calculated affinity spline coefficients a GCRMA pre-processing of Affymetrix Exon microarrays (or other DNA-DNA binding arrays with terminal labeling).
}
\usage{data(affinity.spline.coefs.exon)}
\format{
  The format is:
 Named num [1:15] -0.031403 -0.01374 0.04836 -0.026573 -0.000623 ...
 - attr(*, "names")= chr [1:15] "A" "A" "A" "A" ...
}
\details{
The spline coefficients were calculated by the \code{base.profiles} function from package \code{gcrma} on a data set consisting of 4 human Affymetrix Exon microarrays to which RNA from drosophila was hybridised. The probes on the human array should thus in principle mainly measure NSB.
}
\references{

  Rainer J, Lelong J, Bindreither D, Mantinger C, Ploner C, Geley S,
  Kofler R. (2012) Research resource: transcriptional response to
  glucocorticoids in childhood acute lymphoblastic leukemia. Mol Endocrinol. 2012 Jan;26(1):178-93.

  http://bioinfo.i-med.ac.at
}
\examples{
data(affinity.spline.coefs.exon)

## maybe str(affinity.spline.coefs.exon) ; plot(affinity.spline.coefs.exon) ...
}
\keyword{datasets}
