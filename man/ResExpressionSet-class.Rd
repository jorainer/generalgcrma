\name{ResExpressionSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ResExpressionSet-class}
\alias{residuals,ResExpressionSet-method}
\alias{newResExpressionSet}

\title{Class \code{"ResExpressionSet"}}
\description{
A class extending the \code{ExpressionSet} class from the \code{Biobase}
package with an additional slot \code{residuals} containing aggregated
probe residual information from the RMA model fitting of the individual
probe sets.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ResExpressionSet", assayData, phenoData, featureData, exprs, ...)}.

Objects of the type \code{ExpressionSet} can be casted to
\code{ResExpressionSet} using the \code{newResExpressionSet} function.
}
\section{Slots}{
  \describe{
    \item{\code{residuals}:}{Object of class \code{"matrix"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{annotation}:}{Object of class \code{"character"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}: see
      documentation for \code{ExpressionSet} objects. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}: see
      documentation for \code{ExpressionSet} objects. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly.
}
\section{Methods}{
  \describe{
    \item{residuals}{\code{signature(object = "ResExpressionSet")}: ... }
	 }
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
\keyword{classes}
