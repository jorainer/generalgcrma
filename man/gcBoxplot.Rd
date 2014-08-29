\name{gcBoxplot}
\Rdversion{1.1}
\alias{plotNucleoCont}
\alias{gcBoxplot,AffyBatch-method}
\alias{gcBoxplot}
\title{
Plot GC dependency of probe intensities
}
\description{
Creates boxplots depicting the dependency of probe intensities on a specific nucleotide content (e.g. G and C content) of the probes' sequences.
}
\usage{
gcBoxplot( x, keep.pattern, exclude.pattern, ... )
%doGcBoxplot( x, sequences, nt=c( "C", "G" ), ylim=NULL ) %is deprecated
plotNucleoCont(x, sequences, nt = c("C", "G"), log_trans= TRUE, center=FALSE, ylim=NULL, jahcol=FALSE)
}
\arguments{
  \item{x}{
  The AffyBatch object, a matrix or a numeric vector with probe intensities.}
  \item{keep.pattern}{pattern to restrict the boxplots to probes of
    specific probe sets. If not defined all pm probes will be used.}
  \item{exclude.pattern}{pattern to be used to exclude probes for
     specific probe sets. If not defined all pm probes will be used.}
  \item{...}{these args are passed to plotNucleoCount}
  \item{sequences}{The probe sequences with the same order and length as x if x is not an AffyBath.}
  \item{nt}{the nucleotide(s) for which the plot should be
     generated. Allowed values are A, C, G and T.}
  \item{log_trans}{log2 transform probe intensity values (Default is TRUE)}
   \item{center}{boxplot centered on overall intensity average (mean(x))
   (Default is FALSE)}
  \item{ylim}{allows to manually specify the ylim, or to calculate ylim from the data (default).}
  \item{jahcol}{use fancy colramp to visualise probes nucleotide content
    (Default is FALSE).}
}
\details{
Creates boxplots with the average intensities for probes with the same G
and C (or other nucleotides) content. If \code{x} is a matrix the probe
sequences have to be submitted using the \code{sequences}
parameter. Intensities in the rows of the matrix have to match with the
sequences of the probes. The parameters \code{keep.pattern} and
\code{exclude.pattern} can be used for e.g. custom CDF environments,
like for the Exon microarrays, to generate the plot exclusively for
background probes (e.g. by defining \code{keep.pattern="bg"} if the
probe set names for background probes can be identified by an
\code{bg}). The function assume data in natural scale as input and will
default apply a log2 transform on it. Make sure to set the log_trans
argument to FALSE if you whish to input already log2 transformed data.
}
\value{
creates boxplots.
}
\references{

  Rainer J, Lelong J, Bindreither D, Mantinger C, Ploner C, Geley S,
  Kofler R. (2012) Research resource: transcriptional response to
  glucocorticoids in childhood acute lymphoblastic leukemia. Mol Endocrinol. 2012 Jan;26(1):178-93.

  http://bioinfo.i-med.ac.at
}
\author{
Johannes Rainer, Julien Lelong
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
