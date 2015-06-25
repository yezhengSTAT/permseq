readAllocate=function(object,outfileLoc="./",outputFormat=NULL,ChipThres=NULL,chipFile=NULL,bowtieIndex=NULL,csemDir=NULL,bowtieDir=NULL, vBowtie=NULL,mBowtie=NULL,pBowtie=NULL,fragL=NULL){

  if(!file.exists(outfileLoc)){
    system(paste('mkdir ',outfileLoc,sep=''))
  }
  setwd(outfileLoc)
  # Check parameters  
  if(sum(tolower(outputFormat)==c("tagalign","bed"))==0){
    stop("outputFormat: Not supported format")
  }
  if(length(ChipThres) != 0){
    if(ChipThres<0 | ChipThres>1000)
      stop("ChipThres must be a number between 0 and 1000.")
  }
 if(length(fragL) != 0){
    if(fragL<0)
      stop("fragL (fragment length) must be a positive value.")
  }  
  if(length(vBowtie) != 0){
    if(vBowtie<0 | vBowtie>3)
      stop("vBowtie (-v option in Bowtie) must be 0, 1, 2 or 3.")
  }
  if(length(mBowtie) != 0){
    if(mBowtie%%1!=0 | mBowtie<0)
      stop("mBowtie (-m option in Bowtie) must be a positive integer.")
  }
  if(length(pBowtie) != 0){
    if(pBowtie%%1!=0 | pBowtie<0)
      stop("pBowtie (-p option in Bowtie) must be a positive integer.")
  }
  
 # If object=NULL: Analize ChIP data without using prior files
 if(length(object)==0){
    outfile=NULL
    for(i in 1:length(chipFile)){
      t1=strsplit(chipFile[i],'/')[[1]]
      t1=t1[length(t1)]
      t=strsplit(t1,paste('.',sub(".*\\.","",t1), sep=''))[[1]][1]
      outfile=c(outfile,t)
    }
   # Run Bowtie
    message( "Info: Aligning reads using BOWTIE..." )
    chipAlign = list()
    chipSAM = c()
    for(i in 1:length(chipFile)){
      FileFormat<-substr(chipFile[i], nchar(chipFile[i]), nchar(chipFile[i])+1)
      if(FileFormat=="a") FileFormat<-"f"
      system(paste(bowtieDir,'/bowtie -',FileFormat,' -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ',bowtieIndex, ' ',
      chipFile[i], ' ',outfile[i],'.sam', " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", outfile[i], "_Bowtie_temp.txt",sep=''))

      chipAlign[[chipFile[i]]] = .summary(outfileLoc, paste("priorProcessChIP_", outfile[i], "_Bowtie_temp.txt", sep = ""))
      chipSAM[i] <- paste(outfileLoc,'/',outfile,'.sam',sep='')
    }
   # Run CSEM
    message( "Info: Allocating multi-reads..." )
    for(i in 1:length(chipFile)){
      cmd = paste(csemDir,'/run-csem --sam  -p ',pBowtie,' ',chipSAM[i],' ',fragL,' ',outfile[i],'_permseq',sep='') 
      system(cmd)
    }
  }else{  # If object!=NULL: Analize ChIP data using prior information
    chipSAM=object[['chipSAM']]
    chipFile=object[['chipName']]
    histoneGrpL = object[['histoneGrpL']]
    pBowtie=object[['bowtieInfo']]$pBowtie
    csemDir=object[['csemDir']]
    fragL=object[['fragL']]
    outfile=NULL
    for(i in 1:length(chipFile)){
      t=chipFile[i]
      outfile=c(outfile,t)
    }
    message( "Info: Allocating multi-reads..." )
    for(i in 1:length(chipFile)){
      cmd=paste(csemDir,'/run-csem --sam  -p ',pBowtie,' --prior ',object[['prior']][i],' ',chipSAM[i],' ',fragL,' ',   
      outfile[i],'_permseq',sep='')
      system(cmd)
    }
  }
  
  chipFormat=NULL
  if(!is.null(outputFormat)){
    message( "Info: Converting BAM to other formats..." )
    if(tolower(outputFormat)=='tagalign'){
      chipFormat='_permseq.tagAlign'
      for(i in 1:length(outfile)){
        cmd=paste(csemDir,'/csem-generate-input --tag-align ',outfile[i],'_permseq.bam',' ',outfile[i],'_permseq',sep='')
        system(cmd)
       if(!is.null(ChipThres)) system(paste("awk '$5 >= ",ChipThres," {print $0}' ",outfile[i],'_permseq.tagAlign',">",
       outfile[i],'_permseq',ChipThres,'.tagAlign',sep=""))
     }
    }else if(tolower(outputFormat)=='bed'){
      chipFormat='_permseq.bed'
      for(i in 1:length(outfile)){
        cmd=paste(csemDir,'/csem-generate-input --bed ',outfile[i],'_permseq.bam',' ',outfile[i],'_permseq',sep='')
        system(cmd)
	if(!is.null(ChipThres)) system(paste("awk '$5 >= ",ChipThres," {print $0}' ",outfile[i],'_permseq.bed',">",
        outfile[i],'_permseq',ChipThres,'.bed',sep=""))
      }
    }else{
      print('Not supported format')
    }
  }
  if(!is.null(outputFormat)){
    chipFormat=paste(outfile,chipFormat,sep='')
  }
  # Summarize information
  if(length(object) == 0){
    bowtieInfo=vector('list', 5)
    names(bowtieInfo)=c("bowtieIndex","bowtieDir","vBowtie","mBowtie","pBowtie")
    bowtieInfo[[1]]=bowtieIndex
    bowtieInfo[[2]]=bowtieDir
    bowtieInfo[[3]]=vBowtie
    bowtieInfo[[4]]=mBowtie
    bowtieInfo[[5]]=pBowtie

    return(new(Class = "Prior",
               chipAlign = chipAlign,
               chipSAM = paste(outfileLoc, "/", outfile, ".sam", sep = ""),
               chipAllocate = paste(outfileLoc, "/", outfile,'_permseq.bam',sep=''),
               chipFormat = paste(outfileLoc, "/", chipFormat,sep = ""),
               outfileLoc = outfileLoc,
               chipName = outfile,
               chipNum = length(chipFile),
               fragL = fragL,
               bowtieInfo = bowtieInfo,
               csemDir = csemDir,
               dataNum = length(chipFile)))

  }else{
    outfile_chip <- object@chipSAM
    chipFileName <- gsub(".*/(.*).sam", "\\1", outfile_chip)
    outfile_chip <- paste(chipFileName, ".sam", sep = "")
    dnaseThres <- object@dnaseThres
    dnaseKnots <- object@dnaseKnots
    
    outfile_chipmean <- paste("chip", 1:length(outfile_chip), "_chipmean", sep="")
    
    posLoc_bychr <- vector('list', length(object@dnaseName))
    names(posLoc_bychr) <- object@dnaseName
    posLoc_bychr[[length(object@dnaseName)]] <- object@posLoc_bychr
    #calculate averaged chip read counts according to object (dnase or histone information)
    if(!file.exists((paste(outfileLoc,'/',names(posLoc_bychr)[1],'_',outfile_chipmean[1],sep='')))){
      .chipMeanCounts(object, posLoc_bychr,rep(outfileLoc,length(outfile_chip)),outfile_chip,outfileLoc,outfile_chipmean)
    }
    
    object['chipAllocate'] = paste(outfileLoc, "/", outfile,'_permseq.bam', sep='')
    object['chipFormat'] = paste(outfileLoc, "/", chipFormat, sep = "")
    object['outfileLoc'] = outfileLoc

    return(object)
  }
  message( "Info: ---- Done! ----" )
}
