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


