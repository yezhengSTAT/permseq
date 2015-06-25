priorHistone_init = function(histoneFile = NULL, histoneName = NULL, fragL, AllocThres = 900, chrList = NULL, capping = 0, outfileLoc = "./", chipFile, bowtieIndex, csemDir, bowtieDir, vBowtie = 2, mBowtie = 99, pBowtie = 8, save.files = TRUE){

  if(length(chipFile)>1){
    message( "More than one ChIP-seq data detected: Replicates must be from the same TF (experiment).")
  }
  # Check parameters
  if(fragL<0){
    stop("fragL (fragment length) must be a positive value.")
  }  
  if(vBowtie<0 | vBowtie>3){
    stop("vBowtie (-v option in Bowtie) must be 0, 1, 2 or 3.")
  }
  if(mBowtie%%1!=0 | mBowtie<0){
    stop("mBowtie (-m option in Bowtie) must be a positive integer.")
  }
  if(pBowtie%%1!=0 | pBowtie<0){
    stop("pBowtie (-p option in Bowtie) must be a positive integer.")
  }    
  if(AllocThres<0 | AllocThres>1000){
    stop("AllocThres must be a number between 0 and 1000.")
  }

  #If the output folder does not exist, we create it automatically
  if(!file.exists(outfileLoc)){
    system(paste('mkdir ',outfileLoc,sep=''))
  }
  
  #build the chromosome reference file using bowtie-inspect
  script <- "genRef.pl"
  Fn.Path <- system.file(file.path("Perl", script), package = "permseq")
  bowtieIndexName <-  gsub(".*/(.*)", "\\1", bowtieIndex)
  chrom.ref <- paste(outfileLoc, "/", bowtieIndexName, ".ref", sep = "")
  CMD <- paste(bowtieDir, "/bowtie-inspect -s ", bowtieIndex, " | perl ", Fn.Path, " ", chrom.ref, sep = "")
  system(CMD, intern = TRUE)

  #change into the output folder directory
  setwd(outfileLoc)
  
  # Sumarize Bowtie Input Parameters Information
  bowtieInfo=vector('list', 5)
  names(bowtieInfo)=c("bowtieIndex","bowtieDir","vBowtie","mBowtie","pBowtie")
  bowtieInfo[[1]]=bowtieIndex
  bowtieInfo[[2]]=bowtieDir
  bowtieInfo[[3]]=vBowtie
  bowtieInfo[[4]]=mBowtie
  bowtieInfo[[5]]=pBowtie

  #If the histone name(s) is not given, use the histone directory name
  if(is.null(histoneName)){
    histoneName=1:length(histoneFile)
  }
  link <- vector('list',length(histoneFile))
  names(link) <- histoneName
  histoneAlign <- list()

  # Process each Histone data - align to the reference genome and give the bowtie summary information
  for(i in 1:length(histoneFile)){
    link[[i]] <- .multihistoneProcess(histoneFile[i], fragL, AllocThres, chrList, chrom.ref, capping=0,outfileLoc = paste(outfileLoc,'/',histoneName[i],'/',sep=''), histoneName[i], bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie,save.files)
    
    if((sub(".*\\.", "", histoneFile[i]) == "bam") || (sub(".*\\.", "", histoneFile[i]) == "bed")){ 
      histoneAlign[[histoneName[i]]] <- list()  ## If Histone file is bam or bed: No aligment information
    }else{
      histoneAlign[[histoneName[i]]] <- .summary(paste(outfileLoc, "/",histoneName[i], sep = ""), "priorProcessHistoneBowtie_temp.txt")
    }
    
  }
  
  # Align each ChIP data
  outfile_chip <- NULL
  chipAlign <- list()
  
  #Get the ChIP-seq data name
  for(i in 1:length(chipFile)){
    t1 <- strsplit(chipFile[i],'/')[[1]]
    t1 <- t1[length(t1)]
    t <- strsplit(t1,paste('.',sub(".*\\.","",t1), sep=''))[[1]][1]
    outfile_chip <- c(outfile_chip,t)
  }
  
  setwd(outfileLoc)
  for(i in 1:length(chipFile)){
    FileFormat <- substr(chipFile[i], nchar(chipFile[i]), nchar(chipFile[i])+1) 
    if(FileFormat=="a") FileFormat<-"f" #bowtie allow fastq and fasta file format

    # Bowtie alignment
    system(paste(bowtieDir,'/bowtie -',FileFormat,' -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ', bowtieIndex, ' ',chipFile[i], ' ',outfile_chip[i],'.sam',  " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", outfile_chip[i],"_Bowtie_temp.txt", sep=''))
    
    # Summarize ChIP alignment information 
    chipAlign[[outfile_chip[i]]] <- .summary(outfileLoc, paste("priorProcessChIP_", outfile_chip[i], "_Bowtie_temp.txt", sep = ""))    
  }
  
  outfile_chipmean <- paste(outfile_chip,'_chipmean_temp',sep='')

  # Calcualte averaged chip counts at different groups
  link.temp <- vector('list',length(histoneName))
  names(link.temp) <- histoneName
  for(i in 1:length(histoneName)){
    link.temp[[i]] <- link[[i]][['posLoc_bychr']]
  }
  .chipMeanCounts(prior = NULL, link.temp, rep(outfileLoc, length(outfile_chip)), paste(outfile_chip, '.sam', sep=''), outfileLoc, outfile_chipmean)

  #prepare for the plot to select the best histone to be used as dnase data
  ylim <- vector('list',length(chipFile))
  for(i in 1:length(chipFile))
    ylim[[i]] <- c(0,0)
  reps <- vector('list',length(histoneFile))
  for(i in 1:length(reps)){
    reps[[i]] <- vector('list',length(chipFile))
    for(j in 1:length(chipFile)){
      reps[[i]][[j]] <- read.table(paste(histoneName[i],'_',outfile_chipmean[j],sep=''))
      ylim[[j]][1] <- min(c(ylim[[j]][1],min(unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,]))))
      ylim[[j]][2] <- max(c(ylim[[j]][2],max(unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,]))))  
    }
  }
  # Generate plot used to select one histone data which will be used as "DNase data" in the following step.
  pdf(paste(outfileLoc,'/','marginal_histone_plot.pdf',sep=''))
  for(j in 1:length(outfile_chip))
    for(i in 1:length(reps)){
      .fitPlot(list(reps[[i]][[j]]), link[[i]][['dnaseThres']], link[[i]][['dnaseKnots']],name = outfile_chip[j], xlab = histoneName[i], ylim=ylim[[j]])
    }
  dev.off()
  
  # Summarizing information and create "Prior" class. The information saved in the "Prior" object will be passed to "priorHistone_multi"
  for(i in 1:length(chipFile)){
    link[['chipSAM']] <- c(link[['chipSAM']],paste(outfileLoc,'/',outfile_chip[i],'.sam',sep=''))
    link[['chipUni']] <- c(link[['chipUni']],paste(outfileLoc,'/',outfile_chip[i],'.sam.uni.bed',sep=''))
  }
  link[['histoneName']] <- histoneName
  result <- vector('list',length(histoneName))
  names(result) <- histoneName
  for(i in 1:length(result)){
    result[[i]] <- vector('list',length(outfile_chip))
    names(result[[i]]) <- outfile_chip
    for(j in 1:length(chipFile)){
      result[[i]][[j]] <- unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,])
    }
  }
  link[['meanchip_histone']] <- result
  link[['fragL']] <- fragL
  link[['bowtieInfo']] <- bowtieInfo
  link[['csemDir']] <- csemDir
  link[['chrom.ref']] <- chrom.ref
  system(paste('rm -rf ',outfileLoc,'/*temp',sep=''))

  #create new "Prior" object
  result = new("Prior",
    chipName = outfile_chip,
    chipSAM = link[['chipSAM']],
    chipUni = link[['chipUni']],
    chipNum = length(chipFile),
    histoneName = histoneName,
    fragL = fragL,
    bowtieInfo = bowtieInfo,
    csemDir = csemDir,
    chrom.ref = chrom.ref,
    histoneAlign = histoneAlign,
    chipAlign = chipAlign,
    histoneNum = length(histoneFile),
    outfileLoc = outfileLoc,
    dataNum = length(histoneName)
    )
  for(i in 1:length(histoneFile)){
    result@dnaseHistone[[histoneName[i]]] <- link[[histoneName[i]]]
  }
  return(result)
}

#plot the histone vs ChIP-seq read counts
.fitPlot = function(reps, dnaseThres, dnaseKnots, name, xlab, ylim){
  par(mar = c(5, 5, 4, 2) + 0.1)
  for (i in 1:length(reps)) {
    a = reps[[i]][2, ]/reps[[i]][1, ]
    matplot(dnaseThres, t(a), col = "black", lwd = 1, pch = 20, ylab = "Average ChIP read count", xlab = paste("Histone (", xlab, ") read count", sep = ""), main=name, ylim=ylim)
    for (j in 1:length(dnaseKnots)) abline(v = dnaseKnots[j], col = "blue", lwd = 2)
    }
}
