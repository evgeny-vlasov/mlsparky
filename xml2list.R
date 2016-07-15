
require("XML")
require("RCurl")

fileLoc <- "http://nmrml.org/examples/3/MMBBI_10M12-CE01-1a.nmrML"


ml2sparky <- function(fileLoc){

fullList <- xmlToList(fileLoc)


spectrumDataArray <- fullList$spectrumList$spectrum1D$spectrumDataArray


spectrumDataValue <- spectrumDataArray$text
compressed <- spectrumDataArray$.attrs[1]
elength <- spectrumDataArray$.attrs[2]
bformat <- spectrumDataArray$.attrs[3]



binaryArrayDecode <- function (b64string, what="double", sizeof=4, compression=c("gzip", "bzip2", "xz", "none") ) {
	if (missing(compression)) {
		compression <- "none"
	}
	## Decode. TODO: Check cvParam about the encoding
	raws <- memDecompress(base64Decode(b64string, "raw"), type=compression)
	result <- readBin(raws, n=length(raws)+1, what=what, size=sizeof, endian = "little")
	
}

what <- switch(bformat,
		Complex128 = "double", # that's because complex128 is misleading
		Complex64 = "double",
		Integer32 = "integer",
		Complex32int = "integer",
		"class java.lang.Integer" = "integer",
		Complex64int = "currentlynotsupported")

compression <- ifelse(compressed=="true", "gzip", "none")

data <<- binaryArrayDecode(spectrumDataValue, what=what, compression=compression)

np <<- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$.attrs["numberOfDataPoints"])

nuc <<- as.character(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$acquisitionNucleus["name"])

sf <<- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$irradiationFrequency["value"])

sw <<- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$sweepWidth["value"])

# center ?
# upShift? downShift?
# noiseEst?
# inFolder?
# writeShifts
# cor?  writeNoise?
	


}

## Create a UCSF format spectrum
## outPath - character string; full file path for the newly created spectrum
## np - numeric; the number of points in each dimension
## nuc - character string; the nucleus names for each dimension
## sf - numeric; the spectrometer frequency
## sw - numeric; the sweep width (Hz) in each dimension (only required if 
##	writeShifts is set to FALSE, or for compatibility with Sparky)
## center - numeric; the center (PPM) of the spectrum in each dimension (only 
##	required if writeShifts is set to FALSE)
## upShift - numeric; the upfield chemical shift for the spectrum.
## downShift - numeric; the downfield chemical shift for the spectrum
## noiseEst - numeric; the noise estimate for the output spectrum, to be 
##	included in the output header
## data - numeric vector or matrix; the data for the spectrum as it would be 
##	returned by ucsf2D().  The first data point corresponds to the downfield-
##	most point for the spectrum. The last data point in the matrix corresponds 
##	to the upfield-most point for the spectrum. In other words, if you traverse 
##	the matrix by rows, data points in the matrix start at the bottom-left of 
##	the spectrum and move up and to the right (as the spectrum would normally be 
##	viewed)
## inFolder - list; data and file parameters for input spectrum.  This should 
##	match the	output format of ucsf2D().  If provided, fields from this list
##	will be used as values for the any other arguments that are not provided.
## writeShifts - logical; if TRUE, the upfield and downfield chemical shifts
##	are included in the header for the output file.
## cor - logical; if TRUE, the correction factor applied to the upfield
##	chemical shift, center chemical shift, and sweep width when the UCSF file 
##	was orignally read by rNMR will be negated when the new UCSF file is output.
##	In this case, the up and downfield chemical shifts must be provided.  This 
##	adjustment will not be applied to the upfield and downfield shifts 
##	themselves, but will be used in calculating the sweep width and center.
## writeNoise - logical; if TRUE, the noise estimate for the spectrum will be
##	included in the header for the output file.
## Note: file parameters for 2D spectra should be passed in vector format with 
##	the indirect dimension first