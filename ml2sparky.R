################################################################################
##                                                                            ##
##     Internal functions for reading nmrML files and converting to .ucsf     ##
##                                                                            ##
################################################################################


# fileLoc -- string, location of the .nmrML file. Full path or URL.
# outPath -- string, full path for the new .ucsf file including file name
# as described in the writeUcsf() function. 
# fullList is a list containing all data parsed from .nmrML file
# nucName -- string, nucleus name from .nmrML
# nuc -- string, nucleus name in .ucsf format
# bf -- numeric, the value from .nmrML "irradiationFrequency"
# sf -- numeric, the value from .nmrML "effectiveExcit.ationField"
# sw -- numeric, sweep width in hertz
# ppme -- numeric, upshift in ppm
# ppms -- numeric, downshift in ppm

ml2sparky <- function(fileLoc, outPath){
	
	#loading XML and RCurl packages
		##load XML package
	tryCatch(suppressMessages(library(XML)), error=function(er)
				stop("rNMR requires XML package", call.=FALSE))
	
		##load RCurl package
	tryCatch(suppressMessages(library(RCurl)), error=function(er)
				stop("rNMR requires RCurl package", call.=FALSE))
	
	
	
	# parsing full XML to list
	fullList <- xmlToList(fileLoc)
	
	# getting spectrum data
	spectrumDataArray <- fullList$spectrumList$spectrum1D$spectrumDataArray
	
	
	spectrumDataValue <- spectrumDataArray$text
	compressed <- spectrumDataArray$.attrs[1]
	elength <- spectrumDataArray$.attrs[2]
	bformat <- spectrumDataArray$.attrs[3]
	
	
	# getting spectrum data and decompessing if needed
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
	
	data <- binaryArrayDecode(spectrumDataValue, what=what, compression=compression)
	
	
	# getting values from list
	np <- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$.attrs["numberOfDataPoints"])
	
	nucName <- as.character(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$acquisitionNucleus["name"])
	
	nuc <- switch(nucName, "hydrogen atom" = "H1", 'H2'='2H', 'C13'='13C', 'N15'='15N', 'P31'='31P', 'F19'='19F')
	
	sf <- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$effectiveExcitationField["value"])
	
	bf <- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$irradiationFrequency["value"])
	
	sw <- as.numeric(fullList$acquisition$acquisition1D$acquisitionParameterSet$DirectDimensionParameterSet$sweepWidth["value"])
	
	ppms <- as.numeric(fullList$spectrumList$spectrum1D$xAxis["startValue"])
	
	ppme <- as.numeric(fullList$spectrumList$spectrum1D$xAxis["endValue"])	
	
		
	# writeUcsf is the inner rNMR function which makes new .ucsf file
	writeUcsf(outPath = outPath, np = np, nuc = nuc, sf = sf, sw = sw,upShift = ppme, downShift = ppms, 
			data = data, cor=FALSE, writeShifts=TRUE, writeNoise=FALSE)	
}
