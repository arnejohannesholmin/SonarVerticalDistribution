
#### Functions: ####
# Function for reading sonar per ping output data:
readAllProfosPP <- function(x, ...) {
    if(isTRUE(file.info(x)$isdir)) {
		x <- list.files(x, full.names=TRUE)
	}
	
	# Read the files and rbind:
	#data <- lapply(x, read.table, sep="", header=TRUE)
	data <- lapply(x, data.table::fread, ...)
	data <- data.table::rbindlist(data)
	
	return(data)
}

### # Function for selecting only schools within a mean Sv interval (volume ### backscattering strength in dB):
### subsetSchools <- function(x, lower = -70, upper = -20) {
### 	subset(x, Sv.mean > lower & Sv.mean < upper)
### }

### # Function to exclude schools:
### excludeSchools <- function(x, schoolID, IDcol="Id") {
### 	exclude <- x[[IDcol]] %in% schoolID
### 	x[!exclude, ]
### } 

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

	


mean0 <- function(x) {
	sum(x, na.rm = TRUE) / length(x)
}



#' Convert PROFOS ppe-files to NMDEchosounder or ICESAcoustic files
#' 
#' @param ppeFiles A vector of the PROFOS ppe files, or the path to the directory holding the files.
#' @param acousticFile The path to the acoustic file to write the vertical distribution from the sonar data to.
#' 
#' @import data.table
#' 
#' @export
#' 
profosPPE2LUF20 <- function(
	ppeFiles, 
	acousticFile, 
	
	cruise = "FisherySonar", 
	platform = NA, 
	acocat = 12, 
	freq = 38000, 
	transceiver = 1, 
	logDistance = 0.1, 
	channelThickness = 10, 
	#maxRange = 600, 
	minRange = 10, 
	transducerDepth = 7, 
	depth = "Center.dep", 
	upperChannelDepth = 0, 
	propellerAngleDegrees = 45, 
	checkNAs = TRUE, 
	schoolSvThreshold = NULL, 
	schoolAreaThreshold = NULL, 
	schoolNPingThreshold = NULL, 
	ow = FALSE
) {
    
    # Read the sonar data:
    sonarData <- readAllProfosPP(ppeFiles, na.strings = "N/A")

	#To check if any NA in data before continue
	if(checkNAs) {
		sapply(sonarData, function(x) sum(is.na(x)))
	}
	
    # Set Sv to NA when outside of the valid range:
    if(length(schoolSvThreshold) == 2) {
        outside <- sonarData[, Sv.mean < schoolSvThreshold[1] | Sv.mean > schoolSvThreshold[2]]
        if(any(outside)) {
            message("The schoolSvThreshold (", paste(schoolSvThreshold, collapse = ", "), ") discarded ", sum(outside, na.rm = TRUE), " out of ", length(outside), " pings.")
        }
        sonarData[outside, Sv.mean := NA]
    }
	# And with the school area range:
    if(length(schoolAreaThreshold) == 2) {
        outside <- sonarData[, Area < schoolAreaThreshold[1] | Area > schoolAreaThreshold[2]]
        if(any(outside)) {
            message("The schoolAreaThreshold (", paste(schoolAreaThreshold, collapse = ", "), ") discarded ", sum(outside, na.rm = TRUE), " out of ", length(outside), " pings.")
        }
        sonarData[outside, Sv.mean := NA]
    }
    # And with the number of pings per school:
    if(length(schoolNPingThreshold) == 2) {
        sonarData[, NPingsPerSchool := .N, by = "Id"]
        outside <- sonarData[, NPingsPerSchool < schoolNPingThreshold[1] | NPingsPerSchool > schoolNPingThreshold[2]]
        if(any(outside)) {
            message("The schoolNPingThreshold (", paste(schoolNPingThreshold, collapse = ", "), ") discarded ", sum(outside, na.rm = TRUE), " out of ", length(outside), " pings.")
        }
        sonarData[outside, Sv.mean := NA]
    }
    
	# Add sv:
	sonarData <- addVolumeBackscatteringCoefficient(sonarData)

	# Add DateTime in POSIX format:
	sonarData <- addDateTime(sonarData)

	# Order the schools by time, since this is not a requirement in the PROFOS output files:
	sonarData <- orderSchoolsByTime(sonarData)

	## Write the LUF20:
	writeSonarLUF20(
		sonarData = sonarData, 
		acousticFile = acousticFile, 
		cruise = cruise, 
		platform = platform, 
		acocat = acocat, 
		freq = freq, 
		transceiver = transceiver, 
		logDistance = logDistance, 
		channelThickness = channelThickness, 
		#maxRange = maxRange, 
		minRange = minRange, 
		transducerDepth = transducerDepth, 
		depth = depth, 
		upperChannelDepth = upperChannelDepth, 
		propellerAngleDegrees = propellerAngleDegrees, 
		ow = ow
	)
}



# Function to get the sonar NASC per ping
getSonarNASC <- function(sonarData, transducerDepth = 7, depthLabel = "Center.dep", propellerAngleDegrees = 45, channelThickness = 10) {
    
    # Get the channels:
    maxDepth <- sonarData[, max(get(depthLabel), na.rm = TRUE)]
    channelIntervals <- seq(0, maxDepth, by = channelThickness)
    sonarData[, channelID := findInterval(get(depthLabel), channelIntervals)]
    
    # The total backscactter per ping:
    totalBackscatterBy <- c("DateTime", "channelID")
    sonarData[, totalBackscatter := sum(Sv.lin * Area * channelThickness, na.rm = TRUE), by = totalBackscatterBy]
    sonarDataUnique <- unique(sonarData, by = totalBackscatterBy)
    
    # Get the area in the specific channel:
    sonarDataUnique[, channelArea := getSonarChannnelArea(
        channelID = channelID, 
        tiltAngleDegrees = abs(Trans.tilt[1]), # This can change between pings, and can be given positive or negative					
        transducerDepth = ..transducerDepth, 
        channelThickness = ..channelThickness, 
        propellerAngleDegrees = ..propellerAngleDegrees
    )]
    
    # Caluclate the NASC as total backscatter divvided by channel area:
    sonarDataUnique[, NASC := 4 * pi * 1852^2 * totalBackscatter / channelArea]
    
    # Average per log distance:
    averageNASCBy <- c("logDistanceID", "channelID")
    sonarDataUnique[, NASC := mean(NASC), by = averageNASCBy]
    sonarDataUnique <- unique(sonarDataUnique, by = averageNASCBy)
    
    return(sonarDataUnique)
}









getSonarChannnelArea <- function(channelID, tiltAngleDegrees, transducerDepth, channelThickness, propellerAngleDegrees) {
	
	# The are of a channel (between channelThickness * c(channelID - 1, channelID)) is calculated from the following scheme:
	# Consider the depth d of a beam at horizontal range h from the transducer positioned at (0, d_s). The angle downward relative to the sea surface is a.
	# Then tan(a) = (d - d_s) / h, and h = (d - d_s) / tan(a)
	horizontalRange <- function(depth, transducerDepth, tiltAngleDegrees) {
		tiltAngleRadians <- tiltAngleDegrees * pi/180
		
		## Add a check for school depths smaller than transducer depth:
		#if(any(depth < transducerDepth)) {
		#    stop("transducerDepth (", transducerDepth, ") cannot be larger than the smallest school depth ("#, min(depth, na.rm = TRUE), ").")
		#}
		
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
		depth  = channelThickness * (channelID - 1), 
		transducerDepth  = transducerDepth, 
		tiltAngleDegrees  = tiltAngleDegrees
	)
	horizontalRange2 <- horizontalRange(
		depth  = channelThickness * channelID, 
		transducerDepth  = transducerDepth, 
		tiltAngleDegrees  = tiltAngleDegrees
	)
	
	
	# Get the channnel horizontal area:
	channelAreaFull <- circularRingArea(horizontalRange1, horizontalRange2)
	
	# Discard the propeller water:
	channelAreaFinal <- channelAreaFull * (360 - propellerAngleDegrees) / 360
	
	return(channelAreaFinal)
}





# Function to convert the sonarData to a table with all channels as columns, and write to NMDEchosounder xml file:
writeSonarLUF20 <- function(
	sonarData,  
	acousticFile, 
	cruise = "FisherySonar", 
	platform = NA, 
	
	acocat = 12, 
	freq = 38000, 
	transceiver = 1, 
	logDistance = 0.1, 
	channelThickness = 10, 
	#maxRange = 600, 
	minRange = 10, 
	transducerDepth = 7, 
	depth = "Center.dep", 
	upperChannelDepth = 0, 
	propellerAngleDegrees = 45, 
	
	ow = FALSE
) {
    
    #### echosounder_dataset ####
    # Create first the header table:
    echosounder_dataset <- data.table::data.table(
        platform = platform, 
        cruise = cruise
    )
    
    
    #### distance ####
    # Subset to one row per ping
    print("head(sonarData$DateTime)")
    print(head(sonarData$DateTime))
    sonarDataPerPing <- unique(sonarData, by = "DateTime")
    print("head(sonarDataPerPing$DateTime)")
    print(head(sonarDataPerPing$DateTime))
    
    # Add the time difference and multiply with the speed to get sailed distance:
    timeDiff <- as.numeric(diff(sonarDataPerPing$DateTime), units = "secs")
    timeDiff <- c(timeDiff, utils::tail(timeDiff, 1))
    sonarDataPerPing[, timeDifference := timeDiff]
    sonarDataPerPing[, sailedDistanceIndividual := Ship.speed * timeDifference / 1852]
    sonarDataPerPing[, sailedDistance := cumsum(sailedDistanceIndividual)]
    
    # Define log-distances:
    sonarDataPerPing[, logDistanceID := floor(sailedDistance / logDistance) + 1]
    
    # Add logDistanceID also for the sonarData for use when averaging NASC in each log distance:
    sonarData <- merge(sonarData, sonarDataPerPing[, c("DateTime", "logDistanceID")], by = "DateTime", all.x = TRUE)
    print("head(sonarData$DateTime)___________")
    print(head(sonarData$DateTime))
    
    
    
    
    # Get the log distance start time and position:
    start_time <- sonarDataPerPing[, .(start_time = DateTime[1]), by = "logDistanceID"]$start_time
    lat_start <- sonarDataPerPing[, .(lat_start = Ship.lat[1]), by = "logDistanceID"]$lat_start
    lon_start <- sonarDataPerPing[, .(lon_start = Ship.lon[1]), by = "logDistanceID"]$lon_start
    lastPing <- sonarDataPerPing[, which.max(DateTime)]
    print(head(start_time[-1]))
    stop_time <- start_time[-1]
    stop_time <- c(stop_time, NA)
    print(head(stop_time))
    stop_time[length(stop_time)] <- sonarDataPerPing$DateTime[lastPing]
    print(head(stop_time))
    stop_time <- as.POSIXct(c(start_time[-1], sonarDataPerPing[lastPing, DateTime]), tz = "UTC")
    print(head(stop_time))
    lat_stop <- c(lat_start[-1], sonarDataPerPing[lastPing, Ship.lat])
    lon_stop <- c(lon_start[-1], sonarDataPerPing[lastPing, Ship.lon])
    
    
    
    
    
    
    
    integrator_dist <- sonarDataPerPing[, .(integrator_dist = sum(sailedDistanceIndividual)), by = "logDistanceID"]$integrator_dist
    log_start <- c(0, utils::head(cumsum(integrator_dist), -1))
    
    pel_ch_thickness <- channelThickness
    
    
    start_time <- format(start_time, format = "%Y-%m-%d %H:%M:%OS3")
    stop_time <- format(stop_time, format = "%Y-%m-%d %H:%M:%OS3")
    
    
    
    # Create the distance table:
    distance <-  data.table::data.table(
        log_start = log_start,
        start_time = start_time,
        stop_time = stop_time,
        integrator_dist = integrator_dist,
        pel_ch_thickness = pel_ch_thickness,
        lat_start = lat_start,
        lat_stop = lat_stop,
        lon_start = lon_start,
        lon_stop = lon_stop
    )
    
    
    #### frequency ####
    # Create the frequency table, skipping "threshold", "num_pel_ch", "num_bot_ch", "min_bot_depth", "max_bot_depth", "upper_interpret_depth", "lower_interpret_depth", "lower_integrator_depth" and "quality", "bubble_corr", as these are not used by StoX:
    frequency <-  data.table::data.table(
        log_start = log_start,
        start_time = start_time,
        freq = freq, 
        transceiver = 1, 
        # Hard code this to 0, as it is not needed, and using transducerDepth as before may cause error in RstoxData::StoxAcoustic() when there are data higher than the transducerDepth:
        #upper_integrator_depth = transducerDepth
        upper_integrator_depth = 0
    )
    
    
    #### ch_type ####
    ch_type <-  data.table::data.table(
        log_start = log_start,
        start_time = start_time,
        freq = freq, 
        transceiver = 1, 
        type = "P"
    )
    
    
    #### sa_by_acocat ####
    sa_by_acocat <-  data.table::data.table(
        logDistanceID = unique(sonarDataPerPing$logDistanceID), # Include this for merging with NASC in the sa table below
        log_start = log_start,
        start_time = start_time,
        freq = freq, 
        transceiver = 1, 
        type = "P", 
        acocat =  acocat
    )
    
    
    #### sa ####
    sonarNASC <- getSonarNASC(sonarData, transducerDepth = transducerDepth, depthLabel = "Center.dep", propellerAngleDegrees = propellerAngleDegrees, channelThickness = channelThickness)
        
    
    
    sonarNASCColsToKeep <- c(
        "logDistanceID", 
        "channelID", 
        "NASC"
    )
    sa <- merge(sa_by_acocat, sonarNASC[, sonarNASCColsToKeep, with = FALSE], by = "logDistanceID", all = TRUE)
    
    data.table::setorderv(sa, c("logDistanceID", "channelID"))
    
    sa[, logDistanceID := NULL]

    setnames(sa, old = c('channelID','NASC'), new = c('ch','sa'))
    
    AcousticData <- list(
        list(
            echosounder_dataset = echosounder_dataset, 
            distance = distance, 
            frequency = frequency, 
            ch_type = ch_type, 
            sa_by_acocat = sa_by_acocat, 
            sa = sa
        )
    )
    
    
    
    
    
    
    
	# Write the file:
	if(!file.exists(dirname(acousticFile))) {
		dir.create(dirname(acousticFile))
	}
    
    require(RstoxData)
    # Write to NMDBiotic xml:
    RstoxData:::WriteAcoustic(
        AcousticData, 
        FileNames = acousticFile, 
        namespaces = "http://www.imr.no/formats/nmdechosounder/v1", 
        encoding ="UTF-8", 
        overwrite = ow
    )

	return(AcousticData)
}





