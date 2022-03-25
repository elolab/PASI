# Input: argument "pathways" is a list of pathway structures
# Output: Returns a list of three pathway lists. The first list is a pathway list including numbers of 
#         neighbours, then similarly number of second neighbours and finally occurrences.

ExtractPathwayStatistics = function(pathways){
  
  # Initialize lists to be returned
  neighbors = list()
  secondneighbors = list()
  occurrences = list()
  
  # Detect entrez ids of each node in each pathway
  pathwayentrez = lapply(pathways, function(x){
    entrez = x$nodeinfo$Entrez
    return(entrez)
  })
  
  # Create named vector of occurrences
  alloccurrences = table(unlist(pathwayentrez))
  
  # Investigate pathways one by one
  for(i in 1:length(pathways)){
    
    if(is.list(pathways[[i]]$relationinfo)){ # NOTE: ok, because is.list returns TRUE for data frames
      
      # Pick initial pathway info
      pathway = pathways[[i]]
      nodes = pathway$nodeinfo$Id
      entrez = pathwayentrez[[i]]
      start = pathway$relationinfo$StartId
      end = pathway$relationinfo$EndId
      firstn = NULL
      secondn = NULL
      
      # Detect first and second neighbors for each node
      for(j in 1:length(nodes)){
        node = nodes[j]
        prelationindices = which(start == node)
        crelationindices = which(end == node)
        
        # Detect node ids of child and grand child nodes
        if(length(prelationindices) > 0){
          children = end[prelationindices]
          gcrelationindices = which(start %in% children)
          if(length(gcrelationindices) > 0){ grandchildren = end[gcrelationindices]
          } else grandchildren = NULL
          
        } else{
          children = NULL
          grandchildren = NULL
        }
        
        # Detect node ids of parent and grand parent nodes
        if(length(crelationindices) > 0){
          parents = start[crelationindices]
          gprelationindices = which(end %in% parents)
          if(length(gprelationindices) > 0){ grandparents = start[gprelationindices]
          } else grandparents = NULL
          
        } else{
          parents = NULL
          grandparents = NULL
        }
        
        # Calculate the number of acceptable first and second neighbors
        fneighbors = unique(c(parents, children))
        sneighbors  = unique(c(grandparents, grandchildren))
        sneighbors = setdiff(sneighbors, fneighbors) # one node can't be 1. and 2. neighbor
        nmbrfneighbors = length(fneighbors)
        nmbrsneighbors = length(sneighbors)
        
        # Save the number of first and second neighbors
        firstn = c(firstn, nmbrfneighbors)
        secondn = c(secondn, nmbrsneighbors)
      }
      
      # Save vectors including occurrences, (nmbr of) neighbor and second neighbor
      neighbors[[i]] = firstn
      secondneighbors[[i]] = secondn
      occurrences[[i]] = alloccurrences[entrez]
      
    } else{
      entrez = pathwayentrez[[i]]
      neighbors[[i]] = rep(0, length(entrez))
      secondneighbors[[i]] = rep(0, length(entrez))
      occurrences[[i]] = alloccurrences[entrez]
    }
    
  }
  
  return(list(fneighbors=neighbors, sneighbors=secondneighbors, occurrences=occurrences))
}
