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
    node = x$nodeinfo
    relation = x$relationinfo
    
    # Keep only nodes of allowed types (drop also relations from/to excluded node types)
    dropnodetypeindex = which(!(node$Type %in% allowednodetypes))
    if(length(dropnodetypeindex) > 0){
      if(is.data.frame(relation)){
        ids = node$Id[dropnodetypeindex]
        dropind = which((relation$StartId %in% ids) | (relation$EndId %in% ids))
        relation = relation[setdiff(1:nrow(relation), dropind),,drop=F]
        if(nrow(relation) < 1) relation = NA
      }
      node = node[-dropnodetypeindex,,drop=F]
    }
    
    if(nrow(node) < 1) return(NA)
    
    # Merge nodes with identical Entrez
    duplicatednodeindex = which(duplicated(node$Entrez))
    if(length(duplicatednodeindex) > 0){
      if(is.data.frame(relation)){
        for(i in duplicatednodeindex){
          id = node[i,"Id"]
          replacement = node$Id[match(node[i,"Entrez"],node$Entrez[1:(i-1)])]
          relation[relation$StartId == id,"StartId"] = replacement
          relation[relation$EndId == id,"EndId"] = replacement
        }
      }
      node = node[-duplicatednodeindex,,drop=F]
    }
    
    if(is.data.frame(relation)){
      
      # Keep only relations of allowed types
      droprelationtypeindex = which(!(relation$Type %in% allowedrelationtypes))
      
      # Drop relations from a node to itself
      dropidentityrelations = which(relation$StartId == relation$EndId)
      
      # Merge relations with the same start and end node
      dropduplicatedrelations = which(duplicated(relation[,c("StartId","EndId"),drop=F]))
      
      droprelationindex = unique(c(droprelationtypeindex, dropidentityrelations, dropduplicatedrelations))
      if(length(droprelationindex) > 0) relation = relation[-droprelationindex,,drop=F]
      if(nrow(relation) < 1) relation = NA
    }
    
    pathway = list(nodeinfo=node, relationinfo=relation, pathwayname=x$pathwayname)
    return(pathway)
  })
  names(pathways) = names(paths)
  
  # Drop pathways that lost all their nodes in cleaning (i.e. they are only NA)
  keepind = unlist(lapply(pathways, is.list))
  pathways = pathways[keepind]
  
  return(pathways)
}

