# Input: 'pathtofiles' is a path to the folder containing pathway files (and nothing else). 
# Output: Returns a named list of pathways. Each pathway contains three elements: a data frame containing node info,
#         a data frame containing relation info (or NA if the pathway contains no relations), and a character, which
#         is the pathway name.
ReadFilePathways = function(pathtofiles){
  
  #
  oldwd = getwd()
  setwd(pathtofiles)
  files = dir(getwd())
  files = grep(".txt", files, value=T)
  
  # Read in pathways
  pathways = lapply(files, function(f){
    
    # Read in pathway components
    rawlines = readLines(f)
    splitlines = strsplit(rawlines, split=" ")
    elementnum = lengths(splitlines)
    
    # Check that they are ok
    nodelineindex = which(elementnum[2:length(elementnum)] == 1) + 1
    relationlineindex = which(elementnum[2:length(elementnum)] == 3) + 1
    relations = as.data.frame(splitlines[relationlineindex])
    if(any(!(unlist(relations[1:2,]) %in% rawlines[nodelineindex]))){
      errormessage = paste(c("File pathway",f,"contains relation elements not listed as nodes."), collapse=" ")
      stop(errormessage)
    }
    if(length(nodelineindex) < 1){
      errormessage = paste(c("File pathway",f,"does not contain any nodes. See manual for file format instructions."), collapse=" ")
      stop(errormessage)
    }
    if(any(!(elementnum[2:length(elementnum)] %in% c(1,3)))){
      errormessage = paste(c("File pathway",f,"should contain rows with either 1 (nodes) or 3 (relations) elements. The first line is the only allowed exception."), collapse=" ")
      stop(errormessage)
    }
    if(any(relations[1,] == relations[2,])){
      errormessage = paste(c("File pathway",f,"should not contain relations from a node to itself."), collapse=" ")
      stop(errormessage)
    }
    if(any(!(relations[3,] %in% c("+","-")))){
      errormessage = paste(c("File pathway",f,"should contain only relations with regulation + (activation) or - (inhibition)."), collapse=" ")
      stop(errormessage)
    }
    
    # Identify node lines and initialize node info
    nodeinfo = as.data.frame(matrix(NA, ncol=4, nrow=length(nodelineindex)))
    colnames(nodeinfo) = c("Id","Entrez","Type","Role") # SHOULD I ADD COMPONENT
    
    # Fill node info
    nodeinfo$Id = as.character(1:nrow(nodeinfo))
    nodeinfo$Entrez = as.character(unlist(splitlines[nodelineindex]))
    splitentrez = strsplit(rawlines[nodelineindex], split="_")
    groupindex = which(lengths(splitentrez) > 1)
    nodeinfo$Type = "gene"
    nodeinfo$Type[groupindex] = "group"
    
    # Identify relation lines and initialize relation info
    if(length(relationlineindex) >= 1){
      relationinfo = as.data.frame(matrix(NA, ncol=4, nrow=length(relationindex)))
      colnames(relationinfo) = c("StartId","EndId","Direction","Role")
      
      # Fill relation info
      relationinfo$StartId = nodeinfo[match(as.character(relations[1,]), nodeinfo$Entrez), "Id"]
      relationinfo$EndId = nodeinfo[match(as.character(relations[2,]), nodeinfo$Entrez), "Id"]
      relationinfo$Direction = "activation"
      relationinfo$Direction[which(relations[3,] == "-")] = "inhibition"
    } else relationinfo = NA
    
    pathway = list(nodeinfo=nodeinfo, relationinfo=relationinfo, pathwayname=rawlines[1])
    return(pathway)
  })
  names(pathways) = gsub(".txt", "", files)
  
  setwd(oldwd)
  return(pathways)
}