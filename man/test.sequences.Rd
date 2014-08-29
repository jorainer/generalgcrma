\name{test.sequences}
\Rdversion{1.1}
\alias{test.sequences}
\docType{data}
\title{
Test sequence data set.
}
\description{
A set of probe sequences for 1000 probes.
}
\usage{data(test.sequences)}
%\format{
%  The format is:
% num [1:1000, 1:2] 35 58 38 31 33 25 31 31 25 37 ...
% - attr(*, "dimnames")=List of 2
%  ..$ : chr [1:1000] "4661971" "1718884" "6129182" "4994872" ...
%  ..$ : chr [1:2] "a" "b"
%}
\examples{
data(test.sequences)
head( test.sequences )
## maybe str(test.data) ; plot(test.data) ...
}
\keyword{datasets}
