\name{gcrma.adjust.nomm}
\Rdversion{1.1}
\alias{gcrma.adjust.nomm}
\title{
GCRMA pre-processing of arrays with custom CDF format
}
\description{
\code{gcrma.adjust.nomm} is a wrapper function for functions in the \code{gcrma} package that allows to GCRMA-pre-process Affymetrix arrays with non-typical CDF package and layout (e.g. Affymetrix Exon arrays with custom CDF files). The function performs all steps of the GCRMA preprocessing: optical background correction, NSB background adjustment, quantile normalisation and RMA summarisation.
}
\usage{
gcrma.adjust.nomm(x, optical.correct = TRUE, verbose = TRUE,
affinity.spline.coefs, NCprobe, k = 0.5, GSB.adjust = TRUE, NC.pattern =
"bg", normalize.quantiles = TRUE, normalize.method="quantiles",
normalize.log2=FALSE, summarize = FALSE, summarize.method="rma", cl,
mc.cores, ret.residual.aggr=FALSE, residual.aggr.fun=function( z ){
mean( abs( z ) ) } )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An AffyBatch object
}
  \item{optical.correct}{
If the data should be optical background adjusted.
}
  \item{verbose}{
Print progress information.
}
  \item{affinity.spline.coefs}{
The affinity spline coefficients for the GCRMA NSB adjustment. The package provides pre-calculated coefficients for human Exon arrays, that can be loaded with data( "affinity.spline.coefs.exon" ). Alternatively they can be calculated from the data using the \code{base.profiles.nc} or \code{base.profiles} from the \code{gcrma} package.
}
  \item{NCprobe}{
Indices of non-specific-background control probes within the PM data (not indices of probes in the matrix returned by the functions \code{intensity} or \code{exprs}, but indices of the probes in the intensity matrix returned by the \code{pm} function). The function tries to find NSB control probes within the data using a pattern search on the probe set ids defined in the CDF environment.
}
  \item{k}{
A tuning factor
}
  \item{GSB.adjust}{
If probe effects in specific binding should be adjusted.
}
  \item{NC.pattern}{
Pattern that is used to search for NCprobes in the CDF environment's probe set names.
}
  \item{normalize.quantiles}{
If the data should be quantiles normalised.
}
  \item{normalize.method}{
Which normalization method should be used. At the moment only
\code{quantiles} normalization is implemented.
}
\item{normalize.log2}{
  Wheter the normalization should be performed in log2 or natural scale.
}
  \item{summarize}{
If the probe intensities should be summarised to probe set intensities (using the \code{rma} function/model). If FALSE an AffyBatch will be returned by the function with the background adjusted probe intensities.
}
\item{summarize.method}{
The summarization method/function to be used. Either rma or plm for standard RMA summarization using the median polish algorithm or summarization with the fitPLM function from the affyPLM package.
}
\item{cl}{ Optional paramter. A cluster object created by the \code{makeCluster} function from the \code{snow} package. If submitted a parallelized version of this function is used. }
\item{mc.cores}{ Optional paramter. Number of cores (CPUs) on which the
  background adjustment should be performed. Provides, similar to
  parameter \code{cl}, a significant speed-up of the
  calculation. Requires that the \code{multicore/parallel} package is
  installed. }
\item{ret.residual.aggr}{Only for summarization using the fitPLM
  function (see parameter \code{summarize.method} above). Whether probe
  residual values should be aggregated per probe set. Note that if set
  to \code{TRUE} a \code{ResExpressionSet} will be returned with the
  aggregated residuals available in slot \code{residuals}.}
\item{residual.aggr.fun}{A function that can be used to aggregate
  (combine/summarize) the probe residuals of a probe set.}
}
\details{
The function estimates the affinity spline coefficients using the expression intensities and sequences of the NSB background probes if affinity.spline.coefs is not submitted.
}
\value{
An ExpressionSet with GCRMA preprocessed probe set intensities (if
summarize=TRUE) or an AffyBatch object with GCRMA background adjusted
probe intensity values (in natural scale). If \code{ret.residual.aggr}
is set to \code{TRUE} a \code{ResExpressionSet} object is returned, that
is essentially identical to an \code{ExpressionSet}, with an additional
slot with a matrix (same dimensions than the expression matrix) of the
aggregated probe residuals per probe set.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gcrma.bkg.correct.matrix}}, \code{\link{gcrma.compute.affinities.local.pmonly}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
