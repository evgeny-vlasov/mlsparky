##	the indirect dimension first
writeUcsf <- function(outPath, np, nuc, sf, sw, center, upShift, downShift, 
		noiseEst, data, inFolder, cor=FALSE, writeShifts=FALSE, writeNoise=FALSE){
	
	## Get output path
	if (missing(outPath))
		outPath <- mySave(defaultextension='.ucsf', title='Save spectrum', 
				filetypes=list('ucsf'='UCSF files'))
	if (!length(outPath) || !nzchar(outPath))
		return(invisible())
	
	## Assign arguments if inFolder is provided
	if (!missing(inFolder)){
		if (missing(np))
			np <- inFolder$file.par$matrix_size
		if (missing(nuc))
			nuc <- inFolder$file.par$nucleus
		if (missing(sf))
			sf <- inFolder$file.par$transmitter_MHz
		if (missing(sw))
			sw <- inFolder$file.par$spectrum_width_Hz
		if (missing(center))
			center <- inFolder$file.par$center_ppm
		if (missing(upShift))
			upShift <- inFolder$file.par$upfield_ppm
		if (missing(downShift))
			downShift <- inFolder$file.par$downfield_ppm
		if (missing(noiseEst))
			noiseEst <- inFolder$file.par$noise_est
		if (missing(data))
			data <- inFolder$data
	}else{
		
		## Check for required arguments
		if (missing(np))
			stop('The number of points in the each dimension must be provided')
		if (missing(nuc))
			stop('The nucleus name for each dimension must be provided')
		if (missing(sf))
			stop('The spectrometer frequency for each dimension must be provided')
		if (missing(upShift) || missing(downShift)){
			if (missing(sw))
				stop('The sweep width for each dimension must be provided')
			if (missing(center))
				stop('The center of the spectrum for each dimension must be provided')
			if (writeShifts || cor)
				stop(paste('The up and downfield chemical shifts for each dimension', 
								'must be provided'))
		}
	}
	
	## Negate correction factor originally applied by rNMR
	if (cor){
		cor <- ((downShift - upShift) / (np - 1))	* -(np %% 2 - 1)
		uf <- upShift - cor
		sw <- (downShift - uf) * sf
		center <- downShift - (sw / sf / 2)
	}else{
		
		## Calculate sweep width and center if not provided
		if (missing(sw) || is.null(sw))
			sw <- (downShift - upShift) * sf
		if (missing(center) || is.null(center))
			center <- downShift - ((downShift - upShift) / 2)
	}
	
	## Calculate tile size
	nDim <- length(np)
	tileDim <- np
	if (nDim == 2){
		size <- (tileDim[1] * tileDim[2] * 4) / 1024
		while (size > 32){
			tileDim <- tileDim / 2
			size <- (round(tileDim[1]) * round(tileDim[2]) * 4) / 1024
		}
	}
	tileDim <- round(tileDim)
	
	## Write main sparky header
	if (!file.exists(dirname(outPath)))
		dir.create(dirname(outPath), recursive=TRUE)
	writeCon <- file(outPath, "w+b")
	writeBin('UCSF NMR', writeCon, size=1)
	writeBin(as.integer(0), writeCon, size=1, endian='big')
	writeBin(as.integer(c(nDim, 1, 0, 2)), writeCon, size=1, endian='big')		
	
	## Write out noise estimate ****This differs from UCSF format
	if (writeNoise){
		writeBin('noiseEst', writeCon, size=1, endian='big')
		writeBin(as.numeric(noiseEst), writeCon, size=4, endian='big')
		writeBin(as.integer(rep(0, (180 - 27))), writeCon, size=1, endian='big')
	}else
		writeBin(as.integer(rep(0, (180 - 14))), writeCon, size=1, endian='big')
	
	## Write axis headers
	for (i in 1:nDim){
		writeBin(as.character(nuc[i]), writeCon, size=1)
		writeBin(as.integer(rep(0, (8 - nchar(nuc[i]) - 1))), writeCon, size=1, 
				endian='big')
		writeBin(as.integer(np[i]), writeCon, size=4, endian='big')
		writeBin(as.integer(rep(0, 4)), writeCon, size=1, endian='big')
		writeBin(as.integer(tileDim[i]), writeCon, size=4, endian='big')
		writeBin(as.numeric(sf[i]), writeCon, size=4, endian='big')
		writeBin(as.numeric(sw[i]), writeCon, size=4, endian='big')
		writeBin(as.numeric(center[i]), writeCon, size=4, endian='big')
		
		## Write out upfield and downfield shifts ****This differs from UCSF format
		if (writeShifts){
			writeBin(as.numeric(upShift[i]), writeCon, size=4, endian='big')
			writeBin(as.numeric(downShift[i]), writeCon, size=4, endian='big')
			writeBin(as.integer(rep(0, (128 - 40))), writeCon, size=1, endian='big')
		}else
			writeBin(as.integer(rep(0, (128 - 32))), writeCon, size=1, endian='big')
	}
	
	## Retile and write data out to file
	if (nDim == 1)
		writeBin(as.numeric(data), writeCon, size=4, endian='big')
	else{
		
		## Get data for new tile
		data <- t(data)
		tpc <- ceiling(np[1] / tileDim[1])
		tpr <- ceiling(np[2] / tileDim[2])
		for (i in 1:tpc){
			for (j in 1:tpr){
				rowNum <- (i - 1) * tileDim[1] + 1
				colNum <- (j - 1) * tileDim[2] + 1
				if (j == tpr)
					colOut <- ncol(data) - colNum + 1
				else
					colOut <- tileDim[2]
				if (i == tpc)
					rowOut <- nrow(data) - rowNum + 1
				else
					rowOut <- tileDim[1]
				outData <- data[rowNum:(rowNum + rowOut - 1), 
						colNum:(colNum + colOut - 1)]
				
				## Pad tiles if necessary
				tileRem <- np %% tileDim
				if (all(tileRem != 0) && j == tpr && i == tpc){
					
					## Pad final tile
					if (colOut == 1){
						outData <- c(outData, rep(0, tileDim[1] - length(outData)))
						outData <- cbind(as.numeric(outData), matrix(0, nrow=tileDim[1], 
										ncol=tileDim[2] - 1))
					}else if (rowOut == 1){
						outData <- c(outData, rep(0, tileDim[2] - length(outData)))
						outData <- rbind(as.numeric(outData), matrix(0, 
										nrow=tileDim[1] - 1, ncol=tileDim[2]))
					}else{
						outData <- rbind(outData, matrix(0, nrow=tileDim[1] - nrow(outData), 
										ncol=ncol(outData)))
						outData <- cbind(outData, matrix(0, nrow=tileDim[1], 
										ncol=tileDim[2] - ncol(outData)))
					}
				}else{
					
					## Pad tile in last column
					if (tileRem[2] && j == tpr){
						if (colOut == 1)
							outData <- cbind(as.numeric(outData), matrix(0, nrow=tileDim[1], 
											ncol=tileDim[2] - 1))
						else
							outData <- cbind(outData, matrix(0, nrow=tileDim[1], 
											ncol=tileDim[2] - ncol(outData)))
					}
					## Pad tile in last row
					if (tileRem[1] && i == tpc){
						if (rowOut == 1)
							outData <- rbind(as.numeric(outData), matrix(0, 
											nrow=tileDim[1] - 1, ncol=tileDim[2]))
						else
							outData <- rbind(outData, matrix(0, 
											nrow=tileDim[1] - nrow(outData), ncol=tileDim[2]))
					}
				}						
				
				## Write out new tile
				writeBin(as.numeric(t(outData)), writeCon, size=4, endian='big')
			}
		}
	}
	writeBin('\n', writeCon, size=1)
	close(writeCon)	
	
	return(outPath)
}

