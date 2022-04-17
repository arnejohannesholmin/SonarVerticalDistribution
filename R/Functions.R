
library(stringi)
library(RstoxData)


#### Functions: ####
# Function for reading sonar per ping output data:
readAllProfosPP <- function(x, ...) {
	if(isTRUE(file.info(x)$isdir)) {
		x <- list.files(x, full.names=TRUE)
	}
	
	# Read the files and rbind:
	#data <- lapply(x, read.table, sep="", header=TRUE)
	data <- lapply(x, data.table::fread, ...)
	do.call(rbind, data)
}

# Function for selecting only schools within a mean Sv interval (volume backscattering strength in dB):
subsetSchools <- function(x, lower = -70, upper = -20) {
	subset(x, Sv.mean > lower & Sv.mean < upper)
}

# Function to exclude schools:
excludeSchools <- function(x, schoolID, IDcol="Id") {
	exclude <- x[[IDcol]] %in% schoolID
	x[!exclude, ]
} 

# Function to add volume backsattering coefficient (linear values):
addVolumeBackscatteringCoefficient <- function(x) {
	x$Sv.lin <- 10^(x$Sv.mean/10)
	x
}

# Function to get DateTime:
addDateTime <- function(x, date = "Date", time = "Time", inputformat = "%Y-%m-%d %H:%M:%OS", outputformat = "%Y-%m-%d %H:%M:%OS2") {
	# Convert the date and time of the LUF11:
	isLFU11 <- !grepl(":", x[[time]][1], fixed = TRUE)
	if(isLFU11) {
		x[[date]] <- convertLUF11date(x[[date]])
		x[[time]] <- convertLUF11time(x[[time]])
	}

	# Paste to one date-time string:
	temp <- paste(x[[date]], x[[time]])
	temp <- as.POSIXct(temp, format=inputformat, tz="UTC")
	x$DateTime <- as.POSIXct(temp, format=outputformat, tz="UTC")
	x
}

# Functions to convert LUF11 date and time to ISO format:
convertLUF11time <- function(x) {
	# Add zeros at the start of the time, since the time has been saved as integer, and hours with only one digit are stripped of the zero in the file:
	x <- stringi::stri_pad_left(x, width = 8, pad = "0")
	# Add colons and dot:
	temp <- paste0(
		substr(x, 1, 2), 
		":", 
		substr(x, 3, 4), 
		":", 
		substr(x, 5, 6), 
		".", 
		substr(x, 7, 8)
	)
	temp
}
convertLUF11date <- function(x) {
	# Add dashes:
	temp <- paste0(
		substr(x, 1, 4), 
		"-", 
		substr(x, 5, 6), 
		"-", 
		substr(x, 7, 8)
	)
	temp
}

# Function for ordering the schools by time:
orderSchoolsByTime <- function(x) {
	x[order(x$DateTime), ]
}

# Function to read list user file 11 (LUF11):
readLUF11 <- function(x, headChunk = 20) {
	readLUF11_one <- function(x, headChunk = 20) {
		l <- readLines(x, headChunk)
		skip <- sum(startsWith(l,"%"))
		read.table(x, header = TRUE, sep = "", skip = skip)
	}
	
	if(isTRUE(file.info(x)$isdir)) {
		x <- list.files(x, full.names=TRUE)
	}
	# Read the files and rbind:
	data <- lapply(x, readLUF11_one, headChunk = headChunk)
	do.call(rbind, data)
}	

# Function to get vertical extent of the sonar sampling volume:
#getSonarHeight <- function(tilt = 5, maxRange = 600, minRange = 30, transducerDepth = 6.6, beamWidth = 6) {
#	upperBeamAngleRadians <- pi / 180 * (tilt - beamWidth / 2)
#	lowerBeamAngleRadians <- pi / 180 * (tilt + beamWidth / 2)
#	midBeamAngleRadians <- pi / 180 * tilt
#	upper <- transducerDepth + minRange * sin(upperBeamAngleRadians)
#	lower <- transducerDepth + maxRange * sin(lowerBeamAngleRadians)
#	lowerMid <- transducerDepth + maxRange * sin(midBeamAngleRadians)
#	height <- lower - upper
#	list(
#		upper = upper, 
#		lower = lower, 
#		height = height, 
#		lowerMid = lowerMid
#	)
#}

sonarChannelName <- function(channelID = 1) {
	paste0("Sonar_", echosounderChannelName(channelID))
}
MS70ChannelName <- function(channelID = 1, channelThickness = 5) {
	paste0("sV_Ch", channelID, "_", (channelID - 1) * channelThickness, ".", channelID * channelThickness, "m")
}







setPOSIXlt <- function(x, att.from) {
	attributes(x) <- attributes(att.from)
	#unlist(as.POSIXlt(x))
	x
}

seqIfValid <- function(x, y) {
	if(y >= x) {
		seq(x, y, 1)
	} 
	else {
		NULL
	}
}


mean0 <- function(x) {
	sum(x, na.rm = TRUE) / length(x)
}




# Function to convert the sonarData to a table with all channels as columns, and write to NMDEchosounder xml file:
writeSonarLUF20_old <- function(
	sonarData, 
	LUF20File, 
	cruise = "FisherySonar", 
	platform = NA, 
	acocat = 12, 
	freq = 38000, 
	transceiver = 1, 
	logDuration = 36, 
	channelThickness = 10, 
	#maxRange = 600, 
	minRange = 10, 
	transducerDepth = 7, 
	beamWidth = 6, 
	depth = "Center.dep", 
	upperChannelDepth = 0, 
	propellerAngleDegrees = 45, 
	channelType = "P", 
	xsd = "1", 
	cores = 1
) {
	
	# Get the upper, mid and lower point of the sonar sampling volume:
	numChannels <- ceiling(max(sonarData[[depth]], na.rm = TRUE) / channelThickness)
	channelIDs <- seq_len(numChannels)
	
	# Get the sonar NASC per channel:
	NASCPerPing <- lapply(channelIDs, getSonarNASCPerPing_oneChannel, sonarData = sonarData, transducerDepth = transducerDepth, depth = depth, propellerAngleDegrees = propellerAngleDegrees, channelThickness = channelThickness)
	NASCPerPing <- RstoxData::mergeByIntersect(NASCPerPing, all = TRUE)
	
	# Merge in the positions:
	NASCPerPing <- RstoxData::mergeByIntersect(NASCPerPing, sonarData[, c("DateTime", "Ship.lon", "Ship.lat", "Ship.speed", "Ship.heading", "Trans.tilt")])
	
	# Create a time sequence to interpolate positions to and to average the NASCs and the speed in:
	DateTimeSeq <- seq(min(NASCPerPing$DateTime), max(NASCPerPing$DateTime) + logDuration, by = logDuration)
	# Set the last time of the DateTimeSeq equal to the actual last time:
	DateTimeSeq[length(DateTimeSeq)] <- max(NASCPerPing$DateTime)
	
	# Find first the time intervals of the sonar times:
	timeInterval <- findInterval(NASCPerPing$DateTime, DateTimeSeq)
	
	# Find sequences of unassigned time intervals, and set first time of these sequences to the last time of the DateTimes in the previous interval, and opositely for the last time of each sequence:
	# Locate sequences of unassigned time intervals:
	unassigned <- which(diff(timeInterval) > 1)
	endInterval <- timeInterval[unassigned]
	startInterval <- timeInterval[unassigned + 1]
	
	# Get the end times of each interval which is followed by an empty interval:
	endTime <- unlist(by(
		NASCPerPing$DateTime[timeInterval %in% endInterval], 
		timeInterval[timeInterval %in% endInterval], 
		tail, 
		1
	))
	# Convert to a vector of POSIX:
	endTime <- setPOSIXlt(endTime, NASCPerPing$DateTime)
	
	startTime <- unlist(by(
		NASCPerPing$DateTime[timeInterval %in% startInterval], 
		timeInterval[timeInterval %in% startInterval], 
		head, 
		1
	))
	# Convert to a vector of POSIX:
	startTime <- setPOSIXlt(startTime, NASCPerPing$DateTime)
	
	# Set the times of the endIntervals to the last time in the corresponding pings associatted to the interval. Add one microsecond to the endTime:
	endTime <- endTime + 1e-6
	startTime <- startTime
	DateTimeSeq[endInterval + 1] <- endTime
	DateTimeSeq[startInterval] <- startTime
	
	# Discard the unassigned intervals:
	browser()
	toRemove <- unlist(mapply(seqIfValid, endInterval + 2, startInterval - 1))
	DateTimeSeq <- DateTimeSeq[-toRemove]
	
	# Now we have the time intervals to average the NASC and speed and to get the start and end times and positions in:
	intervalsData <- data.table::data.table(
		DateTime = DateTimeSeq
	)
	# Interpolate the positions onto the time intervals:
	intervalsData[, Longitude := approx(NASCPerPing$DateTime, NASCPerPing$Ship.lon, DateTime)$y]
	intervalsData[, Latitude := approx(NASCPerPing$DateTime, NASCPerPing$Ship.lat, DateTime)$y]
	
	
	# Get start and stop position and time
	LUF20Data <- data.table::data.table(
		intervalsData[-nrow(intervalsData), ], 
		intervalsData[-1, ]
	)
	data.table::setnames(LUF20Data, paste0(rep(c("Start", "Stop"), each = 3), names(LUF20Data)))
	LUF20Data[, timeInterval := seq_len(nrow(.SD))]
	
	# Add the NASC of each channel:
	NASCPerPing[, timeInterval := findInterval(DateTime, DateTimeSeq, rightmost.closed = TRUE)]
	# Average NASC and vessel speed in each interval:
	NASCCols <- names(NASCPerPing)[startsWith(names(NASCPerPing), "sa..ch")]
	num_pel_ch <- length(NASCCols)
	toAverage <- c(
		"Ship.speed", 
		NASCCols
	)
	averageSpeedAndNASC <- NASCPerPing[, lapply(.SD, mean0), by = "timeInterval", .SDcols = toAverage]
	
	# Insert NAs for 0, as per the LUF20 convension:
	for (j in NASCCols) {
		data.table::set(averageSpeedAndNASC, which(averageSpeedAndNASC[[j]] == 0), j, NA)
	}
		
	# Add the average sa and speed to the LUF20Data:
	LUF20Data <- merge(LUF20Data, averageSpeedAndNASC, by = "timeInterval")
	
	
	# Invent log_start as starting from 0. Here speed is in m/s, so we need to convert to knots:
	LUF20Data[, log_start := cumsum(as.numeric(StopDateTime - StartDateTime, units = "hours") * Ship.speed * 3600 / 1852)]
	
	
	integrator_dist <- diff(LUF20Data$log_start)
	integrator_dist <- c(integrator_dist, utils::tail(integrator_dist, 1))
	
	# Add other data:
	LUF20Data[, report_time := format(Sys.time(), tz = "UTC")]
	LUF20Data[, cruise := ..cruise]
	LUF20Data[, platform := ..platform]
	LUF20Data[, integrator_dist := ..integrator_dist]
	LUF20Data[, pel_ch_thickness := ..channelThickness]
	LUF20Data[, num_pel_ch := ..num_pel_ch]
	LUF20Data[, upper_interpret_depth := ..upperChannelDepth]
	LUF20Data[, upper_integrator_depth := ..upperChannelDepth]
	LUF20Data[, acocat := ..acocat]
	LUF20Data[, freq := ..freq]
	LUF20Data[, type := ..channelType]
	LUF20Data[, transceiver := ..transceiver]
	
	
	# Rename variables to the LUF2o names:
	oldNames <- c(
		"StartDateTime", 
		"StopDateTime",
		"StartLongitude",
		"StartLatitude",
		"StopLongitude",
		"StopLatitude"
	)
	LUF20Names <- c(
		"start_time", 
		"stop_time",
		"lon_start",
		"lat_start",
		"lon_stop",
		"lat_stop"
	)
	data.table::setnames(LUF20Data, oldNames, LUF20Names)
	
	# Convert to milliseconds strings:
	LUF20Data[, start_time := format(start_time, format = "%Y-%m-%d %H:%M:%OS3")]
	LUF20Data[, stop_time := format(stop_time, format = "%Y-%m-%d %H:%M:%OS3")]
	
	# Write the file:
	LUF20File <- Rstox::writeAcousticXML(as.data.frame(LUF20Data), LUF20File, xsd = xsd, cores = cores)
	
	return(LUF20Data)
}









profosPP2LUF20 <- function(
	profosDir, 
	LUF20File = file.path(profosDir, "SonarLUF20.xml"), 
	cruise = "FisherySonar", 
	platform = NA, 
	acocat = 12, 
	freq = 38000, 
	transceiver = 1, 
	logDuration = 36, 
	channelThickness = 10, 
	#maxRange = 600, 
	minRange = 10, 
	transducerDepth = 7, 
	beamWidth = 6, 
	depth = "Center.dep", 
	upperChannelDepth = 0, 
	propellerAngleDegrees = 45, 
	channelType = "P", 
	xsd = "1", 
	cores = 1, 
	checkNAs = TRUE, 
	schoolThreshold = c(-70, -20), 
	ow = FALSE
) {
	
	
	#sonarData <- readAllProfosPP(profosDir)
	sonarData <- readAllProfosPP(profosDir, na.strings = "N/A")

	# This was wrong, as we need all pings:
	#sonarData <- sonarData[complete.cases(sonarData), ]

	#To check if any NA in data before continue
	if(checkNAs) {
		sapply(sonarData, function(x) sum(is.na(x)))
	}
	

	# This was also wrong, as we need all pings:
	# Subset by lower and upper mean sv of the schools:
	#sonarData <- subsetSchools(sonarData, lower = schoolThreshold[1], upper = schoolThreshold[2])
	# Rather set Sv to NA when outside of the vvalid range:
	sonarData[Sv.mean < schoolThreshold[1] | Sv.mean > schoolThreshold[2], Sv.mean := NA]

	sonarData <- addVolumeBackscatteringCoefficient(sonarData)

	# Add DateTime in POSIX format:
	sonarData <- addDateTime(sonarData)

	# Order the schools by time, since this is not a requirement in the PROFOS output files:
	sonarData <- orderSchoolsByTime(sonarData)

	# Write the LUF20:
	if(!ow && file.exists(LUF20File)) {
		stop("The LUF20File exists: ", LUF20File, ". Choose a different file path.")
	}
	writeSonarLUF20(
		sonarData = sonarData, 
		LUF20File = LUF20File, 
		cruise = cruise, 
		platform = platform, 
		acocat = acocat, 
		freq = freq, 
		transceiver = transceiver, 
		logDuration = logDuration, 
		channelThickness = channelThickness, 
		#maxRange = maxRange, 
		minRange = minRange, 
		transducerDepth = transducerDepth, 
		beamWidth = beamWidth, 
		depth = depth, 
		upperChannelDepth = upperChannelDepth, 
		propellerAngleDegrees = propellerAngleDegrees, 
		channelType = channelType, 
		xsd = xsd, 
		cores = cores
	)
}














# Function to get the sonar NASC per ping
getSonarNASCPerPing_oneChannel <- function(channelID, sonarData, transducerDepth = 7, depth = "Center.dep", propellerAngleDegrees = 45, channelThickness = 10) {
	
	
	browser()
	# Get the sonar data in the current channel:
	channelUpper <- (channelID - 1) * channelThickness
	channelLower <- channelID * channelThickness
	
	# Here we need to generate 
	sonarDataCopy <- data.table::copy(sonarData)
	inChannel <- (sonarDataCopy[[depth]] >= channelUpper & sonarDataCopy[[depth]] <= channelLower) %in% TRUE
	sonarDataCopy[!inChannel, Sv.lin := NA]
	
	#sonarData <- subset(sonarData, !(sonarData[[depth]] >= channelUpper & sonarData[[depth]] <= channelLower) %in% FALSE)
	
	
	# Get the area weighted sum of the school sv of each ping, and divide by the total sampling area of the sonar (excluding propeller water):
	# R seems to discard decimal seconds in tapply when the grouping varible is POSIX, so we must unclass the time:
	sv_times_area_per_ping <- tapply(sonarDataCopy$Sv.lin * sonarDataCopy$Area, unclass(sonarDataCopy$DateTime), sum)
	
	# Get the area in the specificc channel:
	
	
	# sonarAreaClean <- sonarAreaFull * (1 - propellerAngleDegrees / 360)
	# The tilt is given as either negative or positive (event changing withing one file) for downwards oriented beams, thus the abs():
	tiltAngleDegrees <- abs(unique(sonarDataCopy, by = "DateTime")[["Trans.tilt"]])
	sonarAreaClean <- getSonarChannnelArea(
		channelID = channelID, 
		tiltAngleDegrees = tiltAngleDegrees, # This can change between pings						
		transducerDepth = transducerDepth, 
		channelThickness = channelThickness, 
		propellerAngleDegrees = propellerAngleDegrees
	)
	
	sv_per_ping <- sv_times_area_per_ping / sonarAreaClean
	# Create a data frame including the sv:
	sonarNASCPerPing_oneChannel <- data.table::data.table(
		DateTime = unique(sonarDataCopy$DateTime), 
		sv = sv_per_ping
	)
	# Add the sa:
	sonarNASCPerPing_oneChannel[, sa := sv * channelThickness]
	
	# Convert to NASC:
	sonarNASCPerPing_oneChannel[, paste0("sa..ch.", channelID) :=  4 * pi * 1852^2 * sa]
	
	# Remove the sv and sa, and keep only the NASC (which is sA):
	sonarNASCPerPing_oneChannel[, sv := NULL]
	sonarNASCPerPing_oneChannel[, sa := NULL]
	
	return(sonarNASCPerPing_oneChannel)
}


getSonarChannnelArea <- function(channelID, tiltAngleDegrees, transducerDepth, channelThickness, propellerAngleDegrees) {
	
	# The are of a channel (between channelThickness * c(channelID - 1, channelID)) is calculated from the following scheme:
	# Consider the depth d of a beam at horizontal range h from the transducer positionned at (0, depth d_s). The angle downward relative to the sea surface is a.
	# Then tan(a) = (d - d_s) / h, and h = (d - d_s) / tan(a)
	horizontalRange <- function(depth, transducerDepth, tiltAngleDegrees) {
		tiltAngleRadians <- tiltAngleDegrees * pi/180
		h <- (depth - transducerDepth) / tan(tiltAngleRadians)
		h[h < 0] <- 0
		return(h)
	}
	circularRingArea <- function(range1, range2) {
		pi * (range2^2 - range1^2)
	}
	
	# Get the horizontal ranges for the dephts of the channel:
	depths <- channelThickness * c(channelID - 1, channelID)
	horizontalRange1 <- horizontalRange(
		depth  = depths[1], 
		transducerDepth  = transducerDepth, 
		tiltAngleDegrees  = tiltAngleDegrees
	)
	horizontalRange2 <- horizontalRange(
		depth  = depths[2], 
		transducerDepth  = transducerDepth, 
		tiltAngleDegrees  = tiltAngleDegrees
	)
	
	
	# Get the channnel hoorizontal area:
	channelAreaFull <- circularRingArea(horizontalRange1, horizontalRange2)
	
	# Discard the propeller water:
	channelAreaFinal <- channelAreaFull * (360 - propellerAngleDegrees) / 360
	
	return(channelAreaFinal)
}


# Functions to get log distance IDs:
getLogDistanceInd_one <- function(transect, x) {
	which(x$DateTime >= transect[1] & x$DateTime < transect[2])
}
getLogDistanceInd <- function(transects, x) {
	lapply(transects, getLogDistanceInd_one, x)
}

# Function to compare the depth distribution of echosounder and sonar by the average NASC:
plotNASC_sonar_echosounder <- function(x, scaleSonar = 1) {
	plot(x$echosounder, type="o")
	points(scaleSonar * x$sonar, type="o", col=2)
}

plotMeanNASC_sonar_echosounder <- function(x, scaleSonar = 1, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, channelThickness = 10) {
	sonars <- c("FisherySonar", "MS70", "MS70Phantom")
	for(sonar in sonars) {
		if(length(x[[sonar]])) {
			x[[sonar]] <- x[[sonar]] * scaleSonar
		}
	}
	#x$FisherySonar <- scaleSonar * x$FisherySonar
	#x$MS70 <- scaleSonar * x$MS70
	#x <- gather(x, key = "Instrument", value = "NASC", "Echosounder", "FisherySonar")
	presentInstruments <- setdiff(names(x), c("ChannelID", "TransectID"))
	x <- gather(x, key = "Instrument", value = "NASC", presentInstruments)
	x$Instrument <- factor(x$Instrument, levels = presentInstruments)
	p <- ggplot(data = x) + 
		geom_path(aes(y = ChannelID, x = NASC, colour = Instrument), size=1) + 
		geom_point(aes(y = ChannelID, x = NASC, colour = Instrument), size=3) + 
		coord_cartesian(xlim = xlim, ylim = ylim)
	
	if(!is.na(x$TransectID[1])) {
		p <- p  + ggtitle(paste("Transect", x$TransectID[1]))
	}
	
	if(length(xlab)) {
		p <- p + xlab(xlab)
	}
	if(length(ylab)) {
		p <- p + ylab(ylab)
	}

	getIntervalString <- function(channelID, channelThickness = 5) {
		max <- channelID * channelThickness
		min <- max -  channelThickness
		paste(min, "-", max, "m")
	}
	
	# Reverse the y axis:
	p <- p + scale_y_continuous(
		breaks = x$ChannelID, 
		trans = "reverse", 
		labels = getIntervalString(x$ChannelID, channelThickness = channelThickness)
	)
	
	p
}

# Function to read the Promus10 format:
readPromus10 <- function(x, sonardepth = 7.5, numChannels = 22, freqCompensation = 1.5) {
	# The MS70 sonar shows a peak in the vertical distribution that is some 10-20 meters deeper than the echosounder and fishery sonar. The report nr 10 contains summed sv over time in range and elevation angle bins:
	# Select only lines starting with "2"
	l <- readLines(x)
	rows <- which(startsWith(l, "2"))
	l <- l[rows]
	s <- lapply(l, strsplit, ",", fixed = TRUE)
	s <- lapply(s, unlist)
	s <- lapply(s, as.numeric)
	s <- do.call(rbind, s)
	rangebins <- seq(s[1, 9], s[1, 10])
	phi <- s[,3]
	channelThickness <- s[1, 8]
	range <- channelThickness * rangebins
	# Changed on 2020-08-14:
	#s <- s[, 9 + rangebins]
	s <- s[, 11 + seq_along(rangebins)]
	z <- outer(phi, range, function(phi, range) range * sin((phi - 90) * pi/180)) + sonardepth
	rownames(z) <- phi
	colnames(z) <- range
	
	depthSv <- data.frame(
		depth = c(z), 
		sv <- c(s)
	)
	
	depthsIntervals <- seq(0, numChannels * channelThickness, channelThickness)
	
	
	depthInd <- findInterval(depthSv$depth, depthsIntervals, all.inside = TRUE)
	
	depthSv$depthInd <- depthInd
	
	depthID <- seq_len(length(depthsIntervals) - 1)
	depthMin <- depthsIntervals[-length(depthsIntervals)]
	depthMax <- depthsIntervals[-1]
	sv <- by(depthSv$sv, depthSv$depthInd, mean)
	NASC <- rep(NA, numChannels)
	NASC[as.numeric(names(sv))] <- unclass(sv) * channelThickness * freqCompensation
	#NASC <-  unclass(sv) * channelThickness * freqCompensation
		
	
	data.frame(
		depthID = depthID, 
		depthMin = depthMin, 
		depthMax = depthMax, 
		NASC = NASC)
}














# Function to convert the sonarData to a table with all channels as columns, and write to NMDEchosounder xml file:
writeSonarLUF20 <- function(
	sonarData,  
	LUF20File, 
	cruise = "FisherySonar", 
	platform = NA, 
	acocat = 12, 
	freq = 38000, 
	transceiver = 1, 
	logDuration = 36, 
	channelThickness = 10, 
	#maxRange = 600, 
	minRange = 10, 
	transducerDepth = 7, 
	beamWidth = 6, 
	depth = "Center.dep", 
	upperChannelDepth = 0, 
	propellerAngleDegrees = 45, 
	channelType = "P", 
	xsd = "1", 
	cores = 1
) {
	
	# Get the upper, mid and lower point of the sonar sampling volume:
	numChannels <- ceiling(max(sonarData[[depth]], na.rm = TRUE) / channelThickness)
	channelIDs <- seq_len(numChannels)
	
	# Get the sonar NASC per channel:
	NASCPerPing <- lapply(channelIDs, getSonarNASCPerPing_oneChannel, sonarData = sonarData, transducerDepth = transducerDepth, depth = depth, propellerAngleDegrees = propellerAngleDegrees, channelThickness = channelThickness)
	NASCPerPing <- RstoxData::mergeByIntersect(NASCPerPing, all = TRUE)
	
	# Merge in the positions:
	NASCPerPing <- RstoxData::mergeByIntersect(NASCPerPing, sonarData[, c("DateTime", "Ship.lon", "Ship.lat", "Ship.speed", "Ship.heading", "Trans.tilt")])
	
	browser()
	
	
	# Create a time sequence to interpolate positions to and to average the NASCs and the speed in:
	DateTimeSeq <- seq(min(NASCPerPing$DateTime), max(NASCPerPing$DateTime) + logDuration, by = logDuration)
	# Set the last time of the DateTimeSeq equal to the actual last time:
	DateTimeSeq[length(DateTimeSeq)] <- max(NASCPerPing$DateTime)
	
	# Find first the time intervals of the sonar times:
	timeInterval <- findInterval(NASCPerPing$DateTime, DateTimeSeq)
	

	# Now we have the time intervals to average the NASC and speed and to get the start and end times and positions in:
	intervalsData <- data.table::data.table(
		DateTime = DateTimeSeq
	)
	# Interpolate the positions onto the time intervals:
	intervalsData[, Longitude := approx(NASCPerPing$DateTime, NASCPerPing$Ship.lon, DateTime)$y]
	intervalsData[, Latitude := approx(NASCPerPing$DateTime, NASCPerPing$Ship.lat, DateTime)$y]
	
	
	# Get start and stop position and time
	LUF20Data <- data.table::data.table(
		intervalsData[-nrow(intervalsData), ], 
		intervalsData[-1, ]
	)
	data.table::setnames(LUF20Data, paste0(rep(c("Start", "Stop"), each = 3), names(LUF20Data)))
	LUF20Data[, timeInterval := seq_len(nrow(.SD))]
	
	# Add the NASC of each channel:
	NASCPerPing[, timeInterval := findInterval(DateTime, DateTimeSeq, rightmost.closed = TRUE)]
	# Average NASC and vessel speed in each interval:
	NASCCols <- names(NASCPerPing)[startsWith(names(NASCPerPing), "sa..ch")]
	num_pel_ch <- length(NASCCols)
	toAverage <- c(
		"Ship.speed", 
		NASCCols
	)
	averageSpeedAndNASC <- NASCPerPing[, lapply(.SD, mean0), by = "timeInterval", .SDcols = toAverage]
	
	# Insert NAs for 0, as per the LUF20 convension:
	for (j in NASCCols) {
		data.table::set(averageSpeedAndNASC, which(averageSpeedAndNASC[[j]] == 0), j, NA)
	}
		
	# Add the average sa and speed to the LUF20Data:
	LUF20Data <- merge(LUF20Data, averageSpeedAndNASC, by = "timeInterval")
	
	
	# Invent log_start as starting from 0. Here speed is in m/s, so we need to convert to knots:
	LUF20Data[, log_start := cumsum(as.numeric(StopDateTime - StartDateTime, units = "hours") * Ship.speed * 3600 / 1852)]
	
	
	integrator_dist <- diff(LUF20Data$log_start)
	integrator_dist <- c(integrator_dist, utils::tail(integrator_dist, 1))
	
	# Add other data:
	LUF20Data[, report_time := format(Sys.time(), tz = "UTC")]
	LUF20Data[, cruise := ..cruise]
	LUF20Data[, platform := ..platform]
	LUF20Data[, integrator_dist := ..integrator_dist]
	LUF20Data[, pel_ch_thickness := ..channelThickness]
	LUF20Data[, num_pel_ch := ..num_pel_ch]
	LUF20Data[, upper_interpret_depth := ..upperChannelDepth]
	LUF20Data[, upper_integrator_depth := ..upperChannelDepth]
	LUF20Data[, acocat := ..acocat]
	LUF20Data[, freq := ..freq]
	LUF20Data[, type := ..channelType]
	LUF20Data[, transceiver := ..transceiver]
	
	
	# Rename variables to the LUF2o names:
	oldNames <- c(
		"StartDateTime", 
		"StopDateTime",
		"StartLongitude",
		"StartLatitude",
		"StopLongitude",
		"StopLatitude"
	)
	LUF20Names <- c(
		"start_time", 
		"stop_time",
		"lon_start",
		"lat_start",
		"lon_stop",
		"lat_stop"
	)
	data.table::setnames(LUF20Data, oldNames, LUF20Names)
	
	# Convert to milliseconds strings:
	LUF20Data[, start_time := format(start_time, format = "%Y-%m-%d %H:%M:%OS3")]
	LUF20Data[, stop_time := format(stop_time, format = "%Y-%m-%d %H:%M:%OS3")]
	
	# Write the file:
	if(!file.exists(dirname(LUF20File))) {
		dir.create(dirname(LUF20File))
	}
	LUF20File <- Rstox::writeAcousticXML(as.data.frame(LUF20Data), LUF20File, xsd = xsd, cores = cores)
	
	return(LUF20Data)
}
