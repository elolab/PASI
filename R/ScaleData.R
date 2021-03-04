# Input: "data" is gene level data,
#        "labels" is an integer vector defining sample groups,
#        "noise" is a noise level detected from function DetectNoise,
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# Output: Returns a scaled version of "data" (the order of samples and genes remains the same)

ScaleData = function(data, labels, noise, score){
  
  # Remove genes that are just noise in all groups
  if(!is.na(noise)){
    realsignal = apply(data, 1, function(x){
      groupsignal = unlist(lapply(split(x, labels), median, na.rm=T))
      keep = T
      if(all(groupsignal < noise)) keep = F
      return(keep)
    })
    data = data[realsignal,]
  }
  
  if(score=="activity"){
    
    # Calculate 2 dimensional robust z-scores (takes ages)
    scaleddata = data
    for(i in 1:1000){
      scaleddata = apply(scaleddata, 2, function(x){
        return((x-median(x, na.rm=T))/mad(x, na.rm=T))
      })
      scaleddata = t(scaleddata)
      scaleddata = apply(scaleddata, 2, function(x){
        return((x-median(x, na.rm=T))/mad(x, na.rm=T))
      })
      scaleddata = t(scaleddata)
    }
    
    # Modify Z-scores to interval [-1,1] with logistic function
    #limitedvalues = 2/(1+exp(-2*scaleddata))-1
    scaleddata = 2/(1+exp(-2*scaleddata))-1
  }
  
  if(score=="deregulation"){
    
    samplemedians = apply(data, 2, median, na.rm=T)
    #zeroind = which(data==0, arr.ind=T)
    #zeroindvector = which(data==0)
    greatestmedian = max(samplemedians)
    differences = greatestmedian-samplemedians
    scaleddata = t(t(data)+differences)
    #data[zeroind] = 0
    
    # Calculate normalization parameters for each gene based on expressing control samples
    controls = scaleddata[,labels==0]
    ka = NULL
    kh = NULL
    for(i in 1:nrow(controls)){
      values = as.numeric(controls[i,])
      values = values[!is.na(values)] # Considere only measured controls
      q = quantile(values)
      normrange = 1.5*(q["75%"]-q["25%"])
      values = values[which(values>(q["25%"]-normrange) & values<(q["75%"]+normrange))] # Drop outliers
      ka = c(ka, mean(values))
      kh = c(kh, sd(values))
    }
    parama = -1*ka/kh
    paramb = 1/kh
    
    # Normalize rows using parameters from previous step
    scaleddata = abs(parama+paramb*scaleddata) + 1
    
    # Considere huge outliers
    quantiles = quantile(unlist(scaleddata)) # quantiles of the whole normalized expressed data
    upperlimit = quantiles["75%"] + (1.5*(quantiles["75%"]-quantiles["25%"])) # definition of outlier
    scaleddata[scaleddata>upperlimit] = upperlimit
  }
  
  toreturn = list(scaleddata = scaleddata)
  
  return(toreturn)
}

