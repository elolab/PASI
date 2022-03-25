# Input: "data" is gene level data,
#        "labels" is an integer vector defining sample groups,
#        "noise" is a noise level detected from function DetectNoise,
#        "score" is the score type for final output (character 'activity' or 'deregulation').
# Output: Returns a scaled version of "data" (the order of samples and genes remains the same, though lowly expressed genes are filtered out)
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
  
  # Reduce median control from all samples -> median control has value 0 in all genes
  controlindex = which(labels == 0)
  medians = apply(data[,controlindex], 1, median, na.rm=T)
  centralizeddata = data - medians
  
  # Scale with logistic function to avoid extremes without altering the order (interval [-1,1])
  scaled = t(apply(centralizeddata, 1, function(r){
    k = -log(2/1.5-1) / (quantile(r[controlindex],probs=seq(0,1,by=0.1), na.rm=T)["90%"]) # 90% quantile of controls will get value 75% of maximum 2 (= 1.5)
    temp = 2/(1+exp(-k*r))                                    
    return(temp)
  }))
  scaled = scaled - 1 # logistic function of interval [0,2] - 1 -> interval [-1,1]
  rownames(scaled) = rownames(data)
  colnames(scaled) = colnames(data)
  
  # If deregulation scores are wanted, use absolute values
  if(score == "deregulation") scaled = abs(scaled)
  
  return(scaled)
}

