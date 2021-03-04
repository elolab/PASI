# INPUT: "data" is the gene expression data (rows: Entrez genes, cols: samples),
#        "grouplabels" is an integer vector indicating samples groups (control: 0),
#        "pathwayadress" is path to .xml pathways downloaded from KEGG,
#        "datatype" is should be either "microarray" or "rnaseq",
#        "noisedefault" is the default cutoff (numeric value) for real signal, and
#        "score" defines whether the returned values reflect activity (default) or deregulation.
# OUTPUT: Returns a data frame (cols: samples, rows: pathways) of pathway scores. 

PASI = function(data, grouplabels, pathwayadress=NULL, datatype="microarray", 
                          noisedefault="automatic", score="activity"){
  
  oldwd = getwd()
  
  # Check and initialize input parameters
  parameters = CheckInput(data, grouplabels, pathwayadress, datatype, noisedefault, score)
  originalgenedata = parameters[[1]]
  noisedefault = parameters[[2]]
  
  # Detect noise
  if(is.na(noisedefault)){
    noiselevel = noisedefault
  } else noiselevel = DetectNoise(originalgenedata, noisedefault)
  
  # Detect which KEGG human pathways can be accessed via API (requires internet connection)
  # NOTE: modify the suffix in link for other organisms
  rawapipathways = strsplit(readLines("http://rest.kegg.jp//list/pathway/hsa"), split="\t")
  apipathways = unlist(lapply(rawapipathways, function(x){return(gsub("path:", "", x[1]))}))
  
  # Access KEGG human pathways via API (requires internet connection)
  conprefix = "http://rest.kegg.jp//get/"
  apiconnections = unlist(lapply(apipathways, function(x){paste(conprefix, "/kgml", sep=x)}))
  rawpathwaysautomatic = ReadKEGGPathways(apiconnections)
  
  # Read pathways from file
  if(is.null(pathwayadress)){
    rawpathwaysfile = NULL
  } else{
    setwd(pathwayadress)
    files = dir(pathwayadress)
    files = grep(".xml", files, value=T)
    rawpathwaysfile = ReadKEGGPathways(files)
  }
  
  # Combine file and api pathways
  rawpathways = c(rawpathwaysautomatic, rawpathwaysfile) # possibly duplicated names?
  
  # Determine relations's and nodes's natural states when pathway active
  rawpathways = DetermineRoles(rawpathways, score)
  
  # Process and "clean" the pathways from duplicated and immeasurable relations/nodes 
  mediumpathways = CleanPathways(rawpathways)
  emptypathways = which(!is.list(mediumpathways))
  if(length(emptypathways) > 0) mediumpathways = mediumpathways[-emptypathways]
  
  # Decide which Entrez represents the nodes of type "gene" 
  pathways = OneEntrezPerGeneNode(mediumpathways, originalgenedata)
  
  # Normalize the measurements
  scaleddatares = ScaleData(originalgenedata, grouplabels, noiselevel, score)
  scaleddata = scaleddatares$scaleddata
  
  # Drop genes that don't appear in any pathway
  pathwaygenes = unlist(lapply(pathways, function(x){
    entrez = x$nodeinfo$nodeentrez
    return(unlist(entrez))
  }))
  genedata = scaleddata[rownames(scaleddata) %in% pathwaygenes, ]
  
  # Calculate initial node values (NA for missing values)
  nodevalues = MeasurementsToNodes(pathways, genedata)
  
  # Drop pathways that are not measured for any sample
  droppathways = unlist(lapply(nodevalues, function(x) all(is.na(unlist(x)))))
  dropindex = which(droppathways)
  if(length(dropindex)!=0){
    nodevalues = nodevalues[-dropindex]
    cat("Following pathways are left out from analysis due to too few measured nodes:")+cat("\n")
    cat(names(pathways)[dropindex], sep="\n")
    pathways = pathways[-dropindex]
  } 
  
  # Detect structural info from pathways
  pathwaystatistics = ExtractPathwayStatistics(pathways)
  topologyvalues = CalculateTopologyFactors(pathways, pathwaystatistics)
  
  # Process the node values with feedback and calculate pathway scores
  if(score=="deregulation"){
    nodevalues = ConsidereFeedBack(nodevalues, pathways, pathwaystatistics$occurrences)
    results = CalculatePathwayValues_Deregulation(nodevalues, topologyvalues)
  }

  # Calculate relation values and pathway scores
  if(score=="activity"){
    relationvalues = MeasurementsToRelations(pathways, genedata)
    results = CalculatePathwayValues_Activity(nodevalues, relationvalues, topologyvalues)
  }
  
  # Add row and column names to pathway results
  rownames(results) = names(pathways)
  colnames(results) = colnames(genedata)
  
  # Write the pathway values into a text file
  outputfile = paste("PASI_", ".txt", sep=as.character(Sys.Date()))
  setwd(oldwd)
  write.table(results, file=outputfile, quote=F, sep="\t")
  
  
  return(results)
}