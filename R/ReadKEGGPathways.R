
# Input: 'files' is a named list of kgml files including KEGG pathway info 
# Output: Returns a named list (pathways) of lists. For each pathway, there are three elements in its list:
#         1) nodeinfo is a data frame (rows: nodes, cols: different characteristics), 2) relationinfo is similar to nodeinfo,
#         but contains information about relations within that pathway, and 3) full name of the pathway (a character).
ReadKEGGPathways = function(files, data){
  
  keggpathways = lapply(files, function(x){
    lines = readLines(x)
    
    ##### Node info
    
    # Detect lines including relevant node information
    nodelines = lines[grep("<entry",lines)]
    nodelinelist = strsplit(nodelines, "\"")
    
    # Initialize node info
    nodeinfo = as.data.frame(matrix(NA, ncol=5, nrow=length(nodelinelist)))
    colnames(nodeinfo) = c("Id","Entrez","Component","Type","Role")
    
    # Collect node info
    for(index in 1:length(nodelinelist)){
      
      i = nodelinelist[[index]]
      
      nodeinfo[index,"Id"] = as.character(i[2])
      nodeinfo[index,"Entrez"] = as.character(i[4])
      type = i[6]
      nodeinfo[index,"Type"] = type
      
      # For 'group' nodes, extract components ARE THESE COMPONENTS USED ANYWHERE?
      if(type != "group"){
        nodeinfo[index,"Component"] = 0       
        
      } else {
        lineindex1 = grep(paste(i, collapse="\""), lines)
        lineindex2 = grep("</entry>", lines)
        lineindex2 = lineindex2[lineindex2 > lineindex1][1]
        componentlines = grep("<component",lines[(lineindex1+2):(lineindex2-1)],value=T)
        components = strsplit(componentlines, "\"")
        comp = lapply(components, function(x){return(x[2])})
        nodeinfo[index,"Component"] = paste(unlist(comp), collapse="_")
      }
    }
    
    # Only one Entrez id for nodes of type 'gene' (the one with the highest mean expression in data is selected)
    geneindex = which(nodeinfo$Type == "gene")
    if(length(geneindex) > 0){
      nodeinfo[geneindex,"Entrez"] = unlist(lapply(strsplit(nodeinfo[geneindex,"Entrez"], " "), function(g){
        genes = intersect(gsub("\\D","",g), rownames(data))
        if(length(genes) > 0){
          meangenes = rowMeans(data[genes,,drop=F], na.rm=T)
          toreturn = genes[which.max(meangenes)]
        } else toreturn = g[1]
        return(toreturn)
        })) 
    }
    
    # Group nodes Entrez ids into a format of Entrez1_Entrez2_Entrez3
    groupindex = which(nodeinfo$Type == "group")
    if(length(groupindex) > 0){
      ids = strsplit(nodeinfo[groupindex,"Component"], split="_")
      nodeinfo[groupindex,"Entrez"] = unlist(lapply(ids, function(g){
        allentrez = nodeinfo[nodeinfo$Id %in% g,"Entrez"]
        paste(sort(allentrez),collapse="_")
      })) 
    }
    
    ##### Relation info
    
    # Detect lines starting and ending relations
    startindex = grep("<relation", lines)
    endindex = grep("</relation>", lines)
    
    if(length(startindex) > 0){
      
      # Initialize relation info vectors to be detected
      relationinfo = as.data.frame(matrix(NA, ncol=7, nrow=length(startindex)))
      colnames(relationinfo) = c("StartId","EndId","Type","Name","Value","Direction","Role")
      relationstart = strsplit(lines[startindex], "\"")
      
      for(j in 1:length(startindex)){
        
        # Pick the main relation info from the first relation line
        startline = relationstart[[j]]
        relationinfo[j,"StartId"] = as.character(startline[2])
        relationinfo[j,"EndId"] = as.character(startline[4])
        relationinfo[j,"Type"] = startline[6]
        
        if((endindex[j]-1) > startindex[j]){
          
          # Initialize secundary info and direction from subrelation lines
          sublines = lines[(startindex[j]+1):(endindex[j]-1)]
          subrelations = strsplit(sublines, "\"")
          prenames = NULL
          prevalues = NULL
          interaction = "activation"
          
          # Pick names and values for each subrelation
          for(k in subrelations){
            prenames = c(prenames, k[2])
            prevalues = c(prevalues, k[4])
          }
          
          # Decide direction from subrelation names
          if(any(c("inhibition", "repression") %in% prenames)){
            interaction = "inhibition"
          }
          if(all(c("inhibition", "activation") %in% prenames)){
            interaction = "activation"
          }
          
          # For exception relations (rare) without any secundary information, set artificial info  
        } else{
          prenames = NA
          prevalues = NA
          interaction = "activation"
        }
        
        # Fill in information extracted from subrelations
        relationinfo[j,"Name"] = paste(prenames, collapse="_")
        relationinfo[j,"Value"] = paste(prevalues, collapse="_")
        relationinfo[j,"Direction"] = interaction
      }
      
    } else relationinfo = NA
    
    ##### Pathway name
    
    # Detect pathway name
    nameline = grep("title", lines, value=T)[1]
    pathwayname = unlist(strsplit(nameline, "\""))[2]
    
    return(list(nodeinfo=nodeinfo, relationinfo=relationinfo, pathwayname=pathwayname))
  })
  
  # Use KEGG identifiers as names (example: hsa05130) 
  todropfromname = "(.xml)|(http://rest.kegg.jp//get/)|(/kgml)"
  names(keggpathways) = gsub(todropfromname, "", files)
  return(keggpathways)
}