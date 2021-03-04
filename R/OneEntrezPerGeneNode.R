# INPUT: Argument "pathways" is a list structure including pathway information (nodes other than type gene or
#        group are already filtered out), "unscaleddata" is a data frame of gene level expression data before 
#        normalization. 
# OUTPUT: Returns a list structure similar to "pathways", but each node of type gene has just one associated
#         Entrez id. If there were originally multiple ids, the one with highest average intencity in
#         "unscaleddata" is selected. For non-measured genes the selection is done by random (doesn't matter)

OneEntrezPerGeneNode = function(pathways, unscaleddata){
  toreturn = lapply(pathways, function(x){
    pathway = x
    entrezlist = list()
    oldentrez = pathway[[1]]$nodeentrez
    
    for(i in 1:length(pathway[[1]]$nodetype)){
      
      if(pathway[[1]]$nodetype[i] == "gene"){
        
        currententrez = intersect(oldentrez[[i]], rownames(unscaleddata))
        if(length(currententrez) > 0){
          controlintensities = unscaleddata[currententrez, , drop=FALSE]
          meanintensities = rowMeans(controlintensities, na.rm=TRUE)
          maxindex = which.max(meanintensities)
          newentrez = currententrez[maxindex]

        } else newentrez = oldentrez[[i]][1] #In case of unmeasured node, the node entrez is randomly selected
        
        entrezlist[[i]] = newentrez

      } else{ # node is a group
        entrezlist[[i]] = oldentrez[[i]]
      }
    }
    
    pathway[[1]]$nodeentrez = entrezlist
    return(pathway)
  })
  
  names(toreturn) = names(pathways)
  return(toreturn)
}