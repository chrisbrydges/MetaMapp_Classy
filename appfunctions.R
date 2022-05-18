library(classyfireR)
library(jsonlite)
library(rcdk)
library(shiny)

getChemSimNet <- function (df, smiles,cids, cutoff=0.7) {
  ndf = df
  require(rcdk)

  # This function parses the SMILES
  smiles_list <- sapply(1:nrow(ndf),function(x){
    rcdk::parse.smiles(as.vector(ndf$SMILES)[x])[[1]]
  })
  
  # This function gets the fingerprint. Where it says type="pubchem", that means that it converts the SMILES into binary, which is always 881 digits long
  fps <- t(sapply(1:nrow(ndf), function(x) {
    gc()
    xy <- 0
    xy <- as.character(rcdk::get.fingerprint(smiles_list[[x]],type="pubchem")) 
    xy
  }))
  
  # Next step is to work out where 1s appear in the fingerprint of each compound
  df.bitmat <- do.call(rbind,lapply(fps,function(x) as.integer(strsplit(x,"")[[1]][1:881]))) # 881 because the fingerprints are 881 digits
  df.bitmat.location <- lapply(1:nrow(df.bitmat), function(x) { which(df.bitmat[x,]==1) })
  # Then work out how many 1s appear as well
  fpsmeans <- sapply(df.bitmat.location, function(x){length(x)})
  
  # Calculate Tanimoto similarity for every compound
  s <- matrix(nrow = length(fpsmeans), ncol = length(fpsmeans))
  for(i in 1:length(fpsmeans)){
    for(j in 1:length(fpsmeans)){
      # How many 1s appear for the two compounds
      a <- fpsmeans[i]
      b <- fpsmeans[j]
      # How many 1s appear in common locations across the two compounds
      c <- length(intersect(df.bitmat.location[[i]], df.bitmat.location[[j]]))
      # This is the formula for Tanimoto
      s[i, j] <- c/(a + b - c)
    }
  }
  # Diagonal of the matrix set to zero because we don't need Tanimoto similarities between a compound and itself
  diag(s) <- 0
  dfmax <- cbind(cids,"tmsim",cids[sapply(1:length(cids), function (k) {which.max(s[k,]) })])
  # Set duplicate values (i.e., those below the diagonal of the matrix) to zero as well
  s[lower.tri(s)] <- 0
  # Get rid of compounds with tmsim less than the cutoff
  chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
  chemsimdf <- rbind(chemsimdf, cbind(cids,"tmsim",""), dfmax )
  
  
  ## Duplicated edges removal, These still work 
  chemsimdf <- chemsimdf[-which(chemsimdf[,3]==""),]
  chemsimdf <- rbind(chemsimdf,cbind(chemsimdf[,3], chemsimdf[,2], chemsimdf[,1]))
  chemsimdf <-  chemsimdf[!duplicated( chemsimdf),]
  chemsimdf <- chemsimdf[order(chemsimdf[,3],decreasing = F),]
  chemsimdf <- chemsimdf[order(chemsimdf[,1],decreasing = F),]
  
  pmids_a <- cids
  
  for (i in 1:length(pmids_a)) {
    sind <- c((max(which(chemsimdf[,1]==pmids_a[i])) +1) :nrow(chemsimdf))
    chemsimdf[,3][which(chemsimdf[,3][sind]==pmids_a[i]) + (sind[1]-1) ] <- "XX"
  }
  chemsimdf <- chemsimdf[which(chemsimdf[,3]!="XX"), ]
  #Write the cytoscape network file as an output
  return(chemsimdf)
}



getKEGGRpairs <- function (df,cids, kegg, smiles, cutoff=0.7,krp) {
  keggids = as.vector(df[kegg][[1]])
  cids = as.vector(df[cids][[1]])
  #finds location of keggids indexed on each column of KRPlinks.txt
  krp.1 <- match(krp[,1],keggids)
  krp.2 <- match(krp[,2],keggids)
  #binds them together
  krp.cbind <- cbind (krp.1,krp.2)
  #removes any rows where there isn't a value for both columns (no interaction present between the data)
  krp.net <- subset(krp.cbind, krp.1!="NA" & krp.2!="NA")
  #finds the appropriate pubchemids for each keggid from the input data
  cid.krp.2 <- cids[krp.net[,2]]
  cid.krp.1 <- cids[krp.net[,1]]
  #creates a table of the same interactions but with chemids intead of keggids
  krp.cid.net <- cbind(cid.krp.1,"krp",cid.krp.2)
  chemsim <- getChemSimNet(df,smiles,cids,cutoff)
  krp.cid.net <- rbind(krp.cid.net,chemsim)
  
  return(list(krp = krp.cid.net, chemsim = chemsim))
}




tables = function(df,name, pubchem){
  exportdf = data.frame(Name = df[name][[1]], Pubchem = df[pubchem][[1]])
  InChkeys = vector()
  class = vector()
  superclass = vector()
  for(i in 1:length(exportdf$Pubchem)){
    InCh = jsonlite::fromJSON(url(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",paste0(exportdf$Pubchem[i]),"/property/InChIKey/JSON")))
    Sys.sleep(0.10)
    InChkeys = append(InChkeys,InCh$PropertyTable$Properties[1,2])
    out = tryCatch(
      {
        classyfire = jsonlite::fromJSON(url(paste0("https://cfb.fiehnlab.ucdavis.edu/entities/", paste0(InCh$PropertyTable$Properties[1,2]),".json")))
        if(is.null(classyfire$class$name)){
          class = c(class,'NA')  
        }else{
          class = c(class,classyfire$class$name)
        }
        if(is.null(classyfire$superclass$name)){
          superclass = c(superclass,'NA')  
        }else{
          superclass = c(superclass,classyfire$superclass$name)
        }
        
      },
      error = function(cond){
        print('error')
        
      },
      warning=function(cond){
        classyfire2 = classification(classyfireR::get_classification(InCh$PropertyTable$Properties[1,2]))
        class = c(class, paste0(classyfire2[3,2]))
        if(nrow(classyfire2)>3){
          superclass = c(superclass, paste0(classyfire2[2,2]))
        }else{
          superclass = c(superclass,'NA')
        }
        return(list(class=class,superclass=superclass))
      }
    )
    if(is.list(out)){
      class = out$class
      superclass = out$superclass
    }
  }
  exportdf$class = class
  exportdf$superclass = superclass
  return(exportdf)
}



fcdirection = function(df,pubchemids,pvalue,foldchange){
  cids <- df[pubchemids][[1]]
  
  pvals <- df[pvalue][[1]]
  fc <- df[foldchange][[1]]
  statres <- do.call(rbind,lapply(1:length(pvals),function(z){ if(pvals[z]<0.05) 
  {if (fc[z]>1){return( c(cids[z],pvals[z],round(fc[z],2),"Up")   )}else{  return( c(cids[z],pvals[z],round(1/fc[z],2),"Down")) }     } else {return(c(cids[z],pvals[z],1,"No Change"))} }))
  colnames(statres) <- c("Pubchem","pvalue","foldchange","direction")
  return(statres)
}




getSmiles = function(df,cids){
  output = vector()
  for(i in 1:length(exportdf$Pubchem)){
    Smiles = jsonlite::fromJSON(url(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",paste0(df$cids[i]),"/property/InChIKey/JSON")))
    output = append(output,Smiles$PropertyTable$Properties[1,2])
  }
  return(output)
}



main = function(file1,IndName,IndCID,IndKEGG,IndSmiles,IndPval,IndFC){#under assumption that pubchemids, name, smiles, keggid(not all), fold change, pvalue are given
  df <- read.csv(paste0(file1))
  if( length(which(table(df$IndCID)>1))>0) { stop("Please remove the duplicate PubChem CIDs") }
  if( length(which(table(df$IndKEGG)>1))>0) { stop("Please remove the duplicate KEGG IDs") }
  if( length(table(is.na(df$IndCID)))==2 ) { stop("The data file is missing PubChem CIDs") }
  if( length(table(is.na(df$IndSmiles)))==2 ) { stop("The data file is missing SMILES codes") }
  showNotification("Starting calculations")
  
  progress = shiny::Progress$new()
  on.exit(progress$close())
  n = 4
  
  
  
  
  krp <- read.table("KRPlinks.txt", sep="\t")
  
  progress$inc(1/n, message = paste("Acquiring compound classes"))
  exportdf = tables(df,IndName,IndCID)
  
  progress$inc(1/n, message = paste("Generating network edges"))
  chemsim = getKEGGRpairs(df,IndCID,IndKEGG,IndSmiles,0.7,krp)
  
  progress$inc(1/n, message = paste("Calculating foldchange directions"))
  statres = fcdirection(df,IndCID,IndPval,IndFC)
  progress$inc(1/n, message = paste("Merging data tables"))
  exportdf = merge(x=exportdf,y=statres,by="Pubchem")
  showNotification("Finished calculations")
  return(list(df = exportdf, krp = chemsim[[1]],chemsim = chemsim[[2]]))
}
#artefacts of manually running files
#result = main("metamapp_example (2).csv","ï..CID","KEGG_ID","SMILES","pvalue","foldchange")
#result$df

