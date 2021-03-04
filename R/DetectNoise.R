
# Input: "data" is a numeric matrix, where samples are columns and rows are probes and
#        "noisedefault" is a default noise limit (depends on data type)
# Output: Returns median noise level over samples

DetectNoise = function(data, noisedefault){ 
  
  # Detect the noise level for each sample separately
  noiselevels = apply(data, 2, function(sample){

    # Initial variables
    den = density(sample, na.rm=T)
    logvals = den$x
    densities = den$y
    cutoff = noisedefault
    
    # Define allowed range for the noise level (default +- 15% of the whole log interval)
    lowerlimit = max(cutoff - 0.15*(max(logvals, na.rm=T)-min(logvals, na.rm=T)), min(logvals, na.rm=T))
    upperlimit = min(cutoff + 0.15*(max(logvals, na.rm=T)-min(logvals, na.rm=T)), max(logvals, na.rm=T))
    lowerindex = which.min(abs(logvals-lowerlimit))
    upperindex = which.min(abs(logvals-upperlimit))
    intervaldensities = densities[lowerindex:upperindex]
    intervallogs = logvals[lowerindex:upperindex]
    
    # Detect all canyons in the allowed noise range
    previousvalues = c(-1, intervaldensities[-length(intervaldensities)])
    nextvalues = c(intervaldensities[-1], -1)
    canyonindices = which((intervaldensities < previousvalues) & (intervaldensities < nextvalues))
    
    # If there is a canyon(s), use the lowest one as a cutoff
    if(length(canyonindices) > 0){
      lowestcanyonindex = which.min(intervaldensities[canyonindices])
      canyonindex = canyonindices[lowestcanyonindex]
      cutoff = intervallogs[canyonindex]
    
    # Otherwise search for merging point of two peaks
    } else{
      mergingpointindex = SearchForMergingPoint(intervaldensities, intervallogs)
      if(!is.na(mergingpointindex)) cutoff = intervallogs[mergingpointindex]
    }

    return(cutoff)
  })
  
  # One noise level is used for all the samples. 
  # It might be more sophistigated to define it separatedly for each sample, 
  # but then the code is less robust
  mediancutoff = median(noiselevels) 
  return(mediancutoff)
}





# Helper script for function DetectNoise
# Input: argument "densities" is a density function values of allowed value range and
#        argument "logvals" is numeric vector of allowed range values for which the density values are calculated
# Output: Returns the index of logvalues (= x-axis in density plot) which corresponds to the noise level or NA

SearchForMergingPoint = function(densities, logvals){
  
  # Initialize the return value (and its regulators) as NA (no merging point found)
  mergingpoint = NA
  aftermax = NA
  beforemax = NA
  
  # Detect the highest density in the allowed range (unique peak or end of the interval)
  maxindex = which.max(densities)
  
  # Calculate slopes at each poin (except the ends of the interval)
  y = densities[3:length(densities)] - densities[1:(length(densities)-2)]
  x = logvals[length(logvals)] - logvals[length(logvals)-2] #constant
  slopes = y/x
  
  # If the maximum density is at peak or at the end of the interval, search merging point before maximum
  # If peak is very close to the beginning of the interval, search only after it
  if(maxindex > 10){
    
    # If maximum is peak (not end of inerval), reduce the observed interval to end at the peak
    bslopes = slopes[1:(maxindex-2)]
    
    # Detect where the increase is speading up (slopes more positive)
    bslopedifferences = bslopes[2:length(bslopes)]-bslopes[1:(length(bslopes)-1)]
    slowing = bslopedifferences < 0
    
    # Ignore short slowing/accelerating sections between two opposite long ones
    lengths = rle(slowing)$lengths
    if(length(lengths) > 1){
      small = (lengths < 10)
      longbefore = (c(11,lengths[1:(length(lengths)-1)]) > 10)
      longafter = (c(lengths[2:length(lengths)],11) > 10)
      toignore = which(small & longbefore & longafter)
      if(length(toignore) > 0){
        if(toignore[1] == 1){
          slowing[1:lengths[1]] = !slowing[1:lengths[1]]
          toignore = toignore[-1]
        }
        for(i in toignore){
          startindex = max(sum(lengths[1:(i-1)]) + 1, 1)
          endindex = startindex + lengths[i] - 1
          slowing[startindex:endindex] = !slowing[startindex:endindex]
        }
      }
    }
  
    # Detect where interesting intervals start and end
    intervalstarts = which(slowing[1:(length(slowing)-1)]==T & slowing[2:length(slowing)]==F)+1 # First F
    if(slowing[1] == F) intervalstarts = c(1, intervalstarts)
    intervalends = which(slowing[1:(length(slowing)-1)]==F & slowing[2:length(slowing)]==T) # Last F
    if(slowing[length(slowing)] == F) intervalends = c(intervalends, length(slowing))
    
    if(length(intervalstarts) > 0){ # In case the whole interval is slowing, there's nothing to find
      
      # Pick the longest accelerating interval and use its middlepoint as mergin candidate
      intervallengths = intervalends - intervalstarts + 1
      maxintervalindex = which.max(intervallengths)
      startindex = intervalstarts[maxintervalindex]
      endindex = intervalends[maxintervalindex]
      middleindex = ceiling((endindex+startindex)/2)
      
      # Accept the middlepoint as merging point if slopes before it are negative enough compared to slopes after it
      bearlyslopemedian = median(bslopes[startindex:middleindex])
      blateslopemedian = median(bslopes[(middleindex+1):(endindex+1)])
      if(blateslopemedian/bearlyslopemedian > 1.5){
        beforermax = middleindex + 1 # +1 because the first point doesn't have a slope calculated
      }
    }
  }
  
  # If the maximum density is at peak or at the beginning of the interval, search merging point after maximum
  # If peak is very close to the end of the interval, search only before it
  if(maxindex < length(densities)-10){
    
    # If maximum is peak (not beginning of inerval), reduce the observed interval to start from the peak
    aslopes = slopes[maxindex:length(slopes)]
    
    # Detect where the decrease is slowing down (slopes less negative)
    aslopedifferences = aslopes[2:length(aslopes)]-aslopes[1:(length(aslopes)-1)]
    slowing = aslopedifferences > 0
    
    # Ignore short slowing/accelerating sections between two opposite long ones
    lengths = rle(slowing)$lengths
    if(length(lengths) > 1){
      small = (lengths < 10)
      longbefore = (c(11,lengths[1:(length(lengths)-1)]) > 10)
      longafter = (c(lengths[2:length(lengths)],11) > 10)
      toignore = which(small & longbefore & longafter)
      if(length(toignore) > 0){
        if(toignore[1] == 1){
          slowing[1:lengths[1]] = !slowing[1:lengths[1]]
          toignore = toignore[-1]
        }
        for(i in toignore){
          startindex = max(sum(lengths[1:(i-1)]) + 1, 1)
          endindex = startindex + lengths[i] - 1
          slowing[startindex:endindex] = !slowing[startindex:endindex]
        }
      }
    }
    
    # Detect where interesting intervals start and end
    intervalstarts = which(slowing[1:(length(slowing)-1)]==F & slowing[2:length(slowing)]==T)+1 # First T
    if(slowing[1] == T) intervalstarts = c(1, intervalstarts)
    intervalends = which(slowing[1:(length(slowing)-1)]==T & slowing[2:length(slowing)]==F) # Last T
    if(slowing[length(slowing)] == T) intervalends = c(intervalends, length(slowing))
    
    if(length(intervalstarts) > 0){ # In case the whole interval is accelerating, there's nothing to find
      
      # Pick the longest slowing interval and use its middlepoint as mergin candidate
      intervallengths = intervalends - intervalstarts + 1
      maxintervalindex = which.max(intervallengths)
      startindex = intervalstarts[maxintervalindex]
      endindex = intervalends[maxintervalindex]
      middleindex = ceiling((endindex+startindex)/2)
      
      # Accept the middlepoint as merging point if slopes before it are negative enough compared to slopes after it
      aearlyslopemedian = median(aslopes[startindex:middleindex])
      alateslopemedian = median(aslopes[(middleindex+1):(endindex+1)])
      if(aearlyslopemedian/alateslopemedian > 1.5){
        aftermax = middleindex + maxindex + 1 # +1 because the first point doesn't have a slope calculated
      }
    }
  }
  
  # Decide which one of the merging point candidates to use, if any
  if(!is.na(beforemax)){
    if(!is.na(aftermax)){
      
      # Decide which one of aftermax and beforemax is better (both available)
      if(blateslopemedian/bearlyslopemedian > aearlyslopemedian/alateslopemedian){ mergingpoint = beforemax
      } else mergingpoint = aftermax
      
    } else mergingpoint = beforemax

  } else{
    if(!is.na(aftermax)) mergingpoint = aftermax
  }
  
  return(mergingpoint)
}
