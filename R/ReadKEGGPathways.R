# Input: Vector "files" include names of all KEGG pathway files (.xml) included in the analysis
# Output: Returns a list structure where each element correspond to one pathway. Each pathway element is also
#         a list structure of length 3: node info (list), relation info (list) and full pathway name (string). 

ReadKEGGPathways = function(files){
  keggpathways = lapply(files, function(x){
    lines = readLines(x)
    
    # Detect lines including relevant node information
    nodelines = lines[grep("<entry",lines)]
    nodelinelist = strsplit(nodelines, "\"")
    
    # Initialize node info
    nodeid = NULL
    entreztemp = NULL
    nodecomponent = list()
    nodetype = NULL
    
    # Collect node info
    for(i in nodelinelist){
      type = i[6]
      if(type != "group"){
        nodecomponent[[length(nodecomponent)+1]] = 0       
        
      } else {
        nodecomponent[[length(nodecomponent)+1]] = FindNodeComponents(i, lines)
      }
      entreztemp = c(entreztemp,i[4])
      nodetype = c(nodetype, type)
      nodeid = c(nodeid,i[2])     
    }
    
    # Finalize node info
    entreztemp2 = FixGroupEntrez(entreztemp, nodeid, nodecomponent)
    nodeentrez = FormatNodeEntrez(entreztemp2, nodetype)
    nodeinfo = list(nodeid=nodeid, nodeentrez=nodeentrez, nodecomponent=nodecomponent, nodetype=nodetype)
    
    # Detect lines starting and ending relations
    startindex = grep("<relation", lines)
    endindex = grep("</relation>", lines)
    
    if(length(startindex) > 0){
      
      # Initialize relation info vectors to be detected
      startnode = NULL
      endnode = NULL
      relationtype = NULL
      relationname = NULL
      relationvalue = NULL
      relationstart = strsplit(lines[startindex], "\"")
      direction = NULL
      
      for(j in 1:length(startindex)){
        
        # Pick the main relation info from the first relation line
        startline = relationstart[[j]]
        startnode = c(startnode, startline[2])
        endnode = c(endnode, startline[4])
        relationtype = c(relationtype, startline[6])
        
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
        relationname = c(relationname, paste(prenames, collapse="_"))
        relationvalue = c(relationvalue, paste(prevalues, collapse="_"))
        direction = c(direction, interaction)
      }
      
      # Construct relation info
      relationinfo = list(startnode=startnode, endnode=endnode, relationname=relationname, 
                          relationtype=relationtype, relationvalue=relationvalue, direction=direction)
      
    } else relationinfo = NA
    
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


# Helper script for ReadKEGGPsthways
# Returns a vector of components related to node (if node type is group, there's more than 1 components)

FindNodeComponents = function(line, lines){
  lineindex1 = grep(paste(line, collapse="\""), lines)
  lineindex2 = grep("</entry>", lines)
  lineindex2 = lineindex2[lineindex2 > lineindex1][1]
  componentlines = grep("<component",lines[(lineindex1+2):(lineindex2-1)],value=T)
  components = strsplit(componentlines, "\"")
  comp = lapply(components, function(x){return(x[2])})
  return(unlist(comp))
}

# Helper script for ReadKEGGPsthways
# Returns a vector of Entrez id's related to the node (smallest number first)

FormatNodeEntrez = function(entrez, nodetype){
  allowednodeindex = which(nodetype %in% c("gene","group"))
  splittedentrez = strsplit(as.character(entrez), " ")
  toreturn = lapply(splittedentrez[allowednodeindex], function(x){
    without = gsub("\\D", "", x)
    without = sort(as.numeric(without))
    return(without)
  })
  splittedentrez[allowednodeindex] = toreturn
  return(splittedentrez)
}

# Helper script for ReadKEGGPathways
# Returns a vector of Entrez related to components. If there are multiple Entrez per component, they are 
# combined "entrez1", "entrez2" -> "entrez1 entrez2"

FixGroupEntrez = function(entrez, id, component){
  entreztemp = NULL
  for(i in 1:length(id)){
    if(0 %in% component[[i]]){ entreztemp = c(entreztemp, entrez[i])
    } else{
      componentindex = which(id %in% component[[i]])
      componententrez = entrez[componentindex]
      combinedentrez = paste(componententrez, collapse=" ")
      entreztemp = c(entreztemp, combinedentrez)
    }
  }
  return(entreztemp)
}