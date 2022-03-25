

#
#
#
#
#
#
PreProcessPathways = function(pathwayadress, data, nodemin, score){
  
  # Detect which KEGG human pathways can be accessed via API (requires internet connection)
  # NOTE: modify the suffix in link for other organisms
  rawapipathways = strsplit(readLines("http://rest.kegg.jp//list/pathway/hsa"), split="\t")
  apipathways = unlist(lapply(rawapipathways, function(x){return(gsub("path:", "", x[1]))}))
  
  # Access KEGG human pathways via API (requires internet connection)
  conprefix = "http://rest.kegg.jp//get/"
  apiconnections = unlist(lapply(apipathways, function(x){paste(conprefix, "/kgml", sep=x)}))
  rawpathwaysautomatic = ReadKEGGPathways(apiconnections, data)
  
  # Read pathways from file
  if(is.null(pathwayadress)){
    rawpathwaysfile = NULL
  } else rawpathwaysfile = ReadFilePathways(pathwayadress)

  # Combine file and api pathways
  rawpathways = c(rawpathwaysautomatic, rawpathwaysfile) # possibly duplicated names?
  
  # Which one first: cleaning or roles?
  
  # Determine relations's and nodes's natural states when pathway active
  rawpathways = DetermineRoles(rawpathways, score)
  
  # Process and "clean" the pathways from duplicated and immeasurable relations/nodes 
  pathways = CleanPathways(rawpathways)
  
  # Drop pathways with less than 'nodemin' nodes
  nodenumber = unlist(lapply(pathways, function(p){nrow(p$nodeinfo)}))
  dropindex = which(nodenumber < nodemin)
  if(length(dropindex) > 0){
    ids = names(pathways)[dropindex]
    names = unlist(lapply(pathways[dropindex], function(p){p$pathwayname}))
    pathways = pathways[-dropindex]
    dropped = paste(ids, names, sep=": ")
    print("Following small pathways excluded (see argument 'nodemin'):")
    cat(dropped, sep="\n")
  }
  
  return(pathways)
}