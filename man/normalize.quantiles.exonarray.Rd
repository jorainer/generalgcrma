\name{normalize.quantiles.exonarray}
\Rdversion{1.1}
\alias{normalize.quantiles.exonarray}
\title{
Quantile normalization for Exon microarrays
}
\description{
Performs a between-GeneChip quantile normalization for Affymetrix Exon microarrays or comparable GeneChips.
}
\usage{
normalize.quantiles.exonarray(x, use.log2=FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{
  The AffyBatch object.
}
\item{use.log2}{If the quantile normalisation should be performed in
  log2 scale rather than natural scale. Note that the returned
  normalized intensities are in natural scale.}
  \item{verbose}{
  If additional messages should be printed.
}
}
\details{
Affymetrix Exon microarrays contain probes that target individual Exons
of gene transcripts. Depending on the type of analysis (e.g. summarise
per transcript probes, or per-gene-probes) probes can be part of several
probe sets (i.e. for transcripts of a gene sharing some exons). This
function performs the quantile preprocessing on the intensities of all
unique pm probes (the normalize.quantiles function from the
preprocessCore package would consider the non-unique set of pm probes,
i.e. for transcript probe set definitions the probe sets targeting
transcripts of the same gene share usually a large number of probes).
}
\value{
An AffyBatch object with quantile normalised intensity values (in
natural scale).
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


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
