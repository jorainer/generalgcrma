generalgcrma
============

An R package to perform GCRMA background adjustment for Affymetrix Tiling, Exon and Gene ST microarrays.

Usage
-----

This packages provides functionality to perform a GCRMA pre-processing for the new generation of Affymetrix microarrays (ST GeneChips). To this end the package requires probe sequences of all probes on the microarray as well as a list of probes that measure non-specific signal.
`R` packages providing this data are available at http://dmp.i-med.ac.at/index.php/en/resources/bioinfo-menu-custom-cdf-jumi for some of the newer generation Affymetrix microarrays.

Note that the `affy` package displays a warning message when trying to read the Gene ST microarray and throws an error when reading Exon ST microarrays claiming that the package is not designed to run analyses on these types of microarrays. Actually, using the `generalgcrma` package and the alternative CDF packages above it is possible to analyze them, so, in order to analyze Exon arrays the source code of the `affy` package should be downloaded, the respective code that throws the error message commented out and this modified package installed.

