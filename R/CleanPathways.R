# Input: list "paths" includes pathway elements (also lists).
# Output: Returns a similar structure to the input, but certain type of nodes (and relations including them)  
#         are removed. Also nodes with identical genes and double relations are merged. 
#         Relations from a node to itself are removed as well.

CleanPathways = function(paths){
  
  # All node types are ortholog, enzyme, reaction, gene, group, compound, map and other
  allowednodetypes = c("gene", "group")
  # All relation types are ECrel, PPrel, GErel, PCrel and maplink
  allowedrelationtypes = c("PCrel","GErel","PPrel")    
  
  # Clean and reduce pathways one by one
  pathways = lapply(paths, function(x){
    
    # Extract relevant node and relation info from the pathway
    nodeinf = x$nodeinfo
    relationinf = x$relationinfo
    nodetypes = nodeinf$nodetype
    nodeids = nodeinf$nodeid
    nodeentrez = nodeinf$nodeentrez
    entrez = unlist(lapply(nodeentrez, function(x){
      melted = paste(x, collapse=" ")
      return(melted)
    }))
    if(is.list(relationinf)){
      relationtypes = relationinf$relationtype
      starts = relationinf$startnode
      ends = relationinf$endnode
    }
    
    # Detect indices of nodes with prohibited node type
    dropnodeindex_type = which(!(nodetypes %in% allowednodetypes))
    
    # Detect indices of nodes with identical entrez
    dropnodeindex_duplicate = which(duplicated(entrez))
    
    # Make corresponding changes into relation from/to the double nodes
    if(length(dropnodeindex_duplicate) > 0 & is.list(relationinf)){ 
      replacementindex = match(entrez[dropnodeindex_duplicate],entrez,nomatch=0)
      oldid = nodeids[dropnodeindex_duplicate]
      newid = nodeids[replacementindex]
      for(i in 1:length(oldid)){
        from = oldid[i]
        to = newid[i]
        starts[which(starts==from)] = to
        ends[which(ends==from)] = to
      }
    }
    
    if(is.list(relationinf)){
      
      # Detect indices of relations with prohibited relation type
      droprelationindex_type = which(!(relationtypes %in% allowedrelationtypes))
      
      # Detect indices of relations from or to prohibited node type
      prohibitednodeids = unique(nodeids[dropnodeindex_type])
      droprelationindex_nodetype = which((starts %in% prohibitednodeids)|(ends %in% prohibitednodeids))
      
      # Detect indices of double relations (the first occurrence of the relation is kept)
      droprelationindex_duplicate = which(duplicated(paste(starts, ends, sep=" ")))
      
      # Detect indices of relations from a node to the same node
      droprelationindex_loop = which(starts == ends)
    }
    
    # Form index vectors of nodes and relations that should be kept in the pathway
    dropnodes = unique(c(dropnodeindex_type, dropnodeindex_duplicate))
    keepnodes = setdiff(1:length(nodeids), dropnodes)
    keeprelations = NULL
    if(is.list(relationinf)){
      droprelations =  unique(c(droprelationindex_type, droprelationindex_nodetype, 
                                droprelationindex_duplicate, droprelationindex_loop))
      keeprelations = setdiff(1:length(relationtypes), droprelations)
    }
    
    # If the pathway doesn't include any nodes, return NA, otherwise return reduced pathway
    if(length(keepnodes) > 0){
      if(is.list(relationinf) & (length(keeprelations) > 0)){
        x$relationinfo = list(startnode = starts[keeprelations], 
                              endnode = ends[keeprelations], 
                              relationname = relationinf$relationname[keeprelations], 
                              relationtype = relationinf$relationtype[keeprelations], 
                              relationvalue = relationinf$relationvalue[keeprelations],
                              relationrole = relationinf$relationrole[keeprelations],
                              direction = relationinf$direction[keeprelations])
      } else x$relationinfo = NA
      x$nodeinfo = list(nodeid = nodeinf$nodeid[keepnodes], 
                        nodeentrez = nodeinf$nodeentrez[keepnodes], 
                        nodecomponent = nodeinf$nodecomponent[keepnodes], 
                        nodetype = nodeinf$nodetype[keepnodes],
                        noderole = nodeinf$noderole[keepnodes])
      return(x)
      
    } else return(NA)
  })
  
  names(pathways) = names(paths)
  return(pathways)
}

