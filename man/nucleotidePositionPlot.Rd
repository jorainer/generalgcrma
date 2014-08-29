\name{nucleotidePositionPlot}
\Rdversion{1.1}
\alias{plotNucleoPos}
\alias{nucleotidePositionPlot,AffyBatch-method}
\alias{nucleotidePositionPlot}
\title{
Nucleotide position dependent intensities
}
\description{
Plots the average intensity of all probes with a specific nucleotide per position.
}
\usage{
nucleotidePositionPlot( x, keep.pattern, exclude.pattern ,... )
plotNucleoPos(x, sequences, log_trans = TRUE, center=FALSE)
}
\arguments{
  \item{x}{The AffyBatch object, a matrix or a numeric vector with probe
    intensities.}
  \item{keep.pattern}{pattern to restrict the plot to probes of specific probe sets. If not defined all pm probes will be used.}
  \item{exclude.pattern}{pattern to be used to exclude probes for specific probe sets. If not defined all pm probes will be used.}
  \item{sequences}{The probe sequences with the same order and length as x if x is not an AffyBath.}
  \item{log_trans}{log2 transform probe intensity values (Default is TRUE)}
  \item{center}{mean probe intensity centered on overall intensity average (mean(x))(Default is FALSE).}
  \item{...}{Additional arguments passed to the \code{plotNucleoPos} function.}
}
\details{
Creates a plot of the averaged intensities of probes with a specific nucleotide at each position along the probe sequence.
}
\value{
Creates a plot
}
\references{

  Rainer J, Lelong J, Bindreither D, Mantinger C, Ploner C, Geley S,
  Kofler R. (2012) Research resource: transcriptional response to
  glucocorticoids in childhood acute lymphoblastic leukemia. Mol Endocrinol. 2012 Jan;26(1):178-93.

  http://bioinfo.i-med.ac.at
}
\author{
Julien Lelong, Johannes Rainer
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
