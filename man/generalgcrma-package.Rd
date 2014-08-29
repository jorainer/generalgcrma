\name{generalgcrma-package}
\Rdversion{1.1}
\alias{generalgcrma-package}
\alias{generalgcrma}
\docType{package}
\title{
GCRMA preprocessing of Affymetrix Exon arrays and other non-typical GeneChips
}
\description{
Affymetrix Tiling arrays or Exon arrays usually have no missmatch (MM) probes and have either custom CDF files or no CDF files at all. The \code{gcrma} package however requires that a CDF package and a probe package are available in order to perform the preprocessing. Although gcrma could be applied without taking MM probes into account a straight forward preprocessing with functions from the gcrma package is sometimes not possible. We have therefore developed a package that provides wrapper functions for the main background adjusting functions from the \code{gcrma} package in order to make GCRMA preprocessing and NSB background correction also available for generic input (e.g. performing only probe intensities and their nucleotide sequences).
}
\details{
\tabular{ll}{
Package: \tab generalgcrma\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2012-09-2\cr
License: \tab LGPL\cr
LazyLoad: \tab yes\cr
}

}
\references{

  Rainer J, Lelong J, Bindreither D, Mantinger C, Ploner C, Geley S,
  Kofler R. (2012) Research resource: transcriptional response to
  glucocorticoids in childhood acute lymphoblastic leukemia. Mol Endocrinol. 2012 Jan;26(1):178-93.

  http://bioinfo.i-med.ac.at
}
\author{
Julien Lelong, Johannes Rainer

Maintainer: Johannes Rainer <johannes.rainer@i-med.ac.at>
}
\keyword{ package }
\seealso{
\code{\link[gcrma]{gcrma}}
}
\examples{
}
