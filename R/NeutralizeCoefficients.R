library(MASS)

# Input: 'data' is a matrix of values (rows: genes/pathways/whatever, cols: samples)
#        'times' is a named (col names in data) numeric vector indicating ages/times at samples
#        'controls' is an integer vector indicating which samples in 'data' are controls
#        'donor' is a character vector telling a donor for each sample in 'data'. Can be NA
# Output: Returns a scaeld version of 'data' so that all values are between -1 and 1. The values
#         have been neutralized from effect of donor (if available) and age.
NeutralizeCoefficients = function(data, info, labels, pathwaygenes, mainfeature){
  
  # Preprocess inputs
  pathwaygenes = intersect(pathwaygenes, rownames(data))
  
  # Extract info for numeric features
  info_numeric = info[,apply(info, 2, function(f){any(!is.na(as.numeric(f)))}),drop=F]
  
  # Extract info for categorical features with and without overlap between case and control groups
  controlindex = which(labels == 0)
  caseindex = setdiff(1:length(labels), controlindex)
  info_char = info[,setdiff(colnames(info),colnames(info_numeric)),drop=F]
  info_shared = info_char[,apply(info_char,2,function(g){length(intersect(g[caseindex],g[controlindex]))>1}),drop=F]
  info_split = info_char[,setdiff(colnames(info_char),colnames(info_shared)),drop=F]
  
  # Sort the features in info_split so that binary ones come last
  if(ncol(info_split) > 1){
    
    #
    classes = apply(info, 2, function(f){length(unique(f))})
    info_split = info_split[,c(setdiff(1:length(classes), which(classes==2)), which(classes==2))]
  }
  
  # Extract case and control data for pathway genes
  controldata = data[pathwaygenes,controlindex]
  casedata = data[pathwaygenes,caseindex]
  
  # initialize coefficient data frame
  cathegorynames = unlist(lapply(colnames(info_char), function(g){paste(g,unique(info[,g]),sep="")}))
  cols = c("Intercept_case", "Intercept_control", colnames(info_numeric), cathegorynames)
  coefficients = matrix(0, nrow=length(pathwaygenes), ncol=length(cols))
  colnames(coefficients) = cols
  rownames(coefficients) = pathwaygenes
  
  ### The data contains multiple sample groups
  if(length(unique(labels)) > 1){
    
    controlinfo = cbind(info_numeric[controlindex,,drop=F],cbind(info_shared[controlindex,,drop=F],info_split[controlindex,,drop=F]))
    colnames(controlinfo) = c(colnames(info_numeric),colnames(info_shared),colnames(info_split))
    controlinfo = controlinfo[,which(apply(controlinfo,2,function(x){length(unique(x))})>1),drop=F]
    
    # Record coefficients if the data contains controls and other samples
    for(p in pathwaygenes){
      
      # Construct data frame (controls) for the model fitting
      expressionandinfo = as.data.frame(cbind(as.numeric(t(controldata[p,,drop=F])), controlinfo))
      colnames(expressionandinfo) = c("Expression",colnames(controlinfo))
      if(ncol(info_split) > 0){ # Add prefixes to groups (important for median e.g. donor getting coefficient 0)
        for(f in colnames(info_split)){
          uniquegroups = unique(expressionandinfo[,f])
          medians = NULL
          for(g in uniquegroups) medians = c(medians, median(expressionandinfo[expressionandinfo[,f]==g,1], na.rm=T))
          prefixes = order(order(abs(medians-median(medians))))
          names(prefixes) = uniquegroups
          expressionandinfo[,f] = paste(prefixes[expressionandinfo[,f]], expressionandinfo[,f], sep="_SPAMTAG")
        }
      } 
      
      # Calculate coefficients based on control samples
      coefs_cntrl = rlm(Expression~., data=expressionandinfo)$coefficients
      names(coefs_cntrl) = gsub("\\d*_SPAMTAG", "", names(coefs_cntrl))
      
      # Record those coefficients
      shared = intersect(names(coefs_cntrl), colnames(coefficients))
      coefficients[p,"Intercept_control"] = coefs_cntrl["(Intercept)"]
      coefficients[p,shared] = coefs_cntrl[shared] 
      
      if(ncol(info_split) > 0){
        
        # Remove the effect from case data
        neutraldata = t(casedata[p,,drop=F])
        if(ncol(info_numeric) > 0){
          for(i in intersect(colnames(info_numeric),names(coefs_cntrl))){
            neutraldata = neutraldata - info_numeric[caseindex,i,drop=F]*coefs_cntrl[i]
          }
        }
        if(ncol(info_shared) > 0){
          for(i in colnames(info_shared)){
            categories = paste(i,info_shared[caseindex,i],sep="")
            neutraldata = neutraldata - t(coefficients[p,categories,drop=F]) 
          }
        } 
        
        # Construct data frame (cases) for the model fitting
        expressionandinfo = cbind(neutraldata, info_split[caseindex,,drop=F])
        colnames(expressionandinfo) = c("Expression",colnames(info_split))
        
        # Add prefixes to groups (important for median e.g. donor getting coefficient 0)
        for(f in colnames(info_split)){
          uniquegroups = unique(expressionandinfo[,f])
          medians = NULL
          for(g in uniquegroups) medians = c(medians, median(expressionandinfo[expressionandinfo[,f]==g,1], na.rm=T))
          prefixes = order(order(abs(medians-median(medians))))
          names(prefixes) = uniquegroups
          expressionandinfo[,f] = paste(prefixes[expressionandinfo[,f]], expressionandinfo[,f], sep="_SPAMTAG")
        }
        
        # Calculate intercept and coefficients not overlapping between cases and controls based on case samples
        if((length(mainfeature) == length(labels)) & (all(is.na(mainfeature[controlindex])))){
          expressionandinfo = cbind(expressionandinfo, mainfeature=mainfeature[caseindex])
        }
        coefs_case = rlm(Expression~., data=expressionandinfo)$coefficients
        names(coefs_case) = gsub("\\d*_SPAMTAG", "", names(coefs_case))
        
        # Record coefficient from case samples
        shared = intersect(names(coefs_case), colnames(coefficients))
        coefficients[p,"Intercept_case"] = coefs_case["(Intercept)"]
        if(length(shared) > 0) coefficients[p,shared] = coefs_case[shared]
      }
    }
  } else{
    
    ### If the data contains only one sample group, the situation is more straight forward
    
    # Collect control info
    controlinfo = cbind(info_numeric, cbind(info_shared,info_split))
    colnames(controlinfo) = c(colnames(info_numeric),colnames(info_shared),colnames(info_split))
    if(length(mainfeature) == length(lables)){
      controlinfo = cbind(controlinfo, mainfeature=mainfeature)
    }
    controlinfo = controlinfo[,which(apply(controlinfo,2,function(x){length(unigue(x))})>1)]
    
    # Fit the model and collect coefficients
    for(p in pathwaygenes){
      
      # Calculate coefficients based on control samples
      expressionandinfo = as.data.frame(cbind(t(controldata[p,,drop=F]), controlinfo))
      colnames(expressionandinfo) = c("Expression",colnames(controlinfo))
      coefs_cntrl = rlm(Expression~., data=expressionandinfo)$coefficients
      
      # Record those coefficients
      shared = intersect(names(coefs_cntrl), colnames(coefficients))
      coefficients[p,"Intercept_control"] = coefs_cntrl["(Intercept)"]
      coefficients[p,shared] = coefs_cntrl[shared]
    }
  }
  
  # Remove the effect of given coefficients to be neutralized (NOTE: intercept is not to be neutralized)
  residuals = data[pathwaygenes,]
  if(ncol(info_numeric) > 0){          # neutralize numerical coefficients
    for(n in colnames(info_numeric)) residuals = residuals - coefficients[,n,drop=F] %*% t(info_numeric[,n,drop=F])
  }
  if(ncol(info_numeric) < ncol(info)){ # neutralize categorical coefficients
    for(m in setdiff(colnames(info),colnames(info_numeric))){
      residuals = residuals - coefficients[,paste(m,info[,m],sep="")]
    }
  }
  data[pathwaygenes,] = residuals
  
  return(data)
}