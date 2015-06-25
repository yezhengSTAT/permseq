priorProcess=function(dnaseFile=NULL, histoneFile=NULL, dnaseName='dnase', histoneName=NULL, fragL, AllocThres = 900, chrList = NULL, capping = 0, outfileLoc = "./", outfile, bowtieIndex, csemDir, bowtieDir, vBowtie = 2, mBowtie = 99, pBowtie = 8, save.files = TRUE){
  
  # Check parameters
  if(fragL<0)
    {
      stop("fragL (fragment length) must be a positive value.")
    }  
  if(vBowtie<0 | vBowtie>3)
    {
      stop("vBowtie (-v option in Bowtie) must be 0, 1, 2 or 3.")
    }
  if(mBowtie%%1!=0 | mBowtie<0)
    {
      stop("mBowtie (-m option in Bowtie) must be a positive integer.")
    }
  if(pBowtie%%1!=0 | pBowtie<0)
    {
      stop("pBowtie (-p option in Bowtie) must be a positive integer.")
    }    
  if(AllocThres<0 | AllocThres>1000)
    {
      stop("AllocThres must be a number between 0 and 1000.")
    }
  if(length(dnaseFile)>1)
    {
      stop("Only one DNase file is acceptable. DNase file needs to be merged prior running the priorProcess function.")
    }   
  if(is.null(histoneName))
    {
      histoneName=1:length(histoneFile)
    }
  if(is.null(dnaseName)) 
    {
      dnaseName='dnase'
    }
  
  #build the reference chromosome file using bowtie-inspect
  if(!file.exists(outfileLoc)) system(paste('mkdir ',outfileLoc,sep=''))
  script <- "genRef.pl"
  Fn.Path <- system.file(file.path("Perl", script), package = "permseq")
  bowtieIndexName <-  gsub(".*/(.*)", "\\1", bowtieIndex)
  chrom.ref <- paste(outfileLoc, "/", bowtieIndexName, ".ref", sep = "")
  CMD <- paste(bowtieDir, "/bowtie-inspect -s ", bowtieIndex, " | perl ", Fn.Path, " ", chrom.ref, sep = "")
  system(CMD, intern = TRUE)
  
  # Sumarize Bowtie Information
  bowtieInfo=vector('list', 5)
  names(bowtieInfo)=c("bowtieIndex","bowtieDir","vBowtie","mBowtie","pBowtie")
  bowtieInfo[[1]]=bowtieIndex
  bowtieInfo[[2]]=bowtieDir
  bowtieInfo[[3]]=vBowtie
  bowtieInfo[[4]]=mBowtie
  bowtieInfo[[5]]=pBowtie

  # if only one DNase data
  if(is.null(histoneFile)&!is.null(dnaseFile)){ 
      # Process DNase data
      dnaseObject=.dnaseProcess(dnaseFile,fragL,AllocThres,chrList,chrom.ref,capping=0,
        outfileLoc=paste(outfileLoc,'/',dnaseName,'/',sep=''),paste(dnaseName,'_',outfile,sep=''),  
        bowtieIndex,csemDir,bowtieDir,vBowtie=2,mBowtie=99,pBowtie=8,save.files)
      if((sub(".*\\.", "", dnaseFile) == "bam") || (sub(".*\\.", "", dnaseFile) == "bed")){
        dnaseAlign = list()  ## If DNase file is bam or bed: No aligment information
      }else{
        dnaseAlign = .summary(paste(outfileLoc,'/',dnaseName, sep = ""), "priorProcessDNaseBowtie_temp.txt")
      }
      return(new("Prior",
                 dnaseKnots = dnaseObject[['dnaseKnots']],                 
                 dnaseThres = dnaseObject[['dnaseThres']],                 
                 posLoc_bychr= dnaseObject[['posLoc_bychr']],                 
                 dnaseAlign = dnaseAlign,                 
                 dataNum = 1,                 
                 dnaseName = dnaseName,                 
                 fragL = fragL,                 
                 bowtieInfo = bowtieInfo,                 
                 csemDir = csemDir,          
                 outfileLoc = outfileLoc,             
                 chrom.ref = chrom.ref))
      
    # If DNase and Histone data    
    }else if(!is.null(histoneFile)&!is.null(dnaseFile)){
      histoneNum = length(histoneFile)
      # Process DNase data
      dnaseObject=.dnaseProcess(dnaseFile,fragL,AllocThres,chrList,chrom.ref,capping=0,
        outfileLoc=paste(outfileLoc,'/',dnaseName,'/',sep=''), paste(dnaseName,'_',outfile,sep=''), 
        bowtieIndex,csemDir,bowtieDir,vBowtie=2,mBowtie=99,pBowtie=8,save.files)

      chrList = dnaseObject[['chrList']]

      if((sub(".*\\.", "", dnaseFile) == "bam") || (sub(".*\\.", "", dnaseFile) == "bed")){
        dnaseAlign = list()  ## If DNase file is bam or bed: No aligment information
      }else{
        dnaseAlign = .summary(paste(outfileLoc,'/',dnaseName, sep = ""), "priorProcessDNaseBowtie_temp.txt")
      }
       
      link=vector('list',length(histoneFile))
      names(link)=histoneName
      histoneAlign = list()
      # Process each Histone data     
      for(i in 1:length(histoneFile))
        {
          link[[i]]=.histoneProcess(histoneFile[i],fragL,AllocThres,chrList,chrom.ref,capping=0, outfileLoc=paste (outfileLoc,'/',histoneName[i],'/',sep=''), paste(histoneName[i],'_',outfile,sep=''),bowtieIndex,csemDir,bowtieDir,vBowtie=2,mBowtie=99,pBowtie=8,save.files)
          if((sub(".*\\.", "", histoneFile[i]) == "bam") || (sub(".*\\.", "", histoneFile[i]) == "bed")){
            histoneAlign[[histoneName[i]]] = list()   ## If Histone file is bam or bed: No aligment information
          }else{
            histoneAlign[[histoneName[i]]] = .summary(paste(outfileLoc, "/",histoneName[i], sep = ""), "priorProcessHistoneBowtie_temp.txt")
          }
        }
      
      #generate multipatterns
      s=scan(chrom.ref,what='char',sep='\n')
      chr_all=strsplit(s[3],' ')[[1]]
      chr_length=as.numeric(strsplit(s[2],' ')[[1]])
      message( "Info: Partitioning genome based on multiple data sets..." )
      for(i in 1:length(chrList)) 
        {
          chr=chrList[i]
          file_input=dnaseObject[['posLoc_bychr']][[chr]]    
          for(h in 1:length(histoneFile)) 
            {
              file_input=paste(file_input, ' ',link[[h]][[chr]],sep='')
            }
          # partition genome according to dnase (dnase and histone data)
          script='multi_data_pattern.pl'
          Fn.Path <- system.file( file.path("Perl",script), package="permseq")
          
          CMD<-paste('perl ',Fn.Path,' ',outfileLoc,'/',chr,'_',outfile,'_positions_cluster.txt',' ', chr_length[which(chr_all==chr)], ' ',file_input,sep='')
          system(CMD)
        }
      link.t=vector('list',length(chrList))
      names(link.t)=chrList
      for(i in chrList) 
        {
          link.t[[i]]=paste(outfileLoc,"/",i,"_",outfile, "_positions_cluster.txt",sep="")
        }
      
      return(new("Prior",
                 dnaseName = dnaseName,                 
                 dnaseKnots = dnaseObject[['dnaseKnots']],             
                 dnaseThres = dnaseObject[['dnaseThres']],                 
                 posLoc_bychr= link.t,             
                 dataNum = 1 + length(histoneFile),                 
                 dnaseAlign = dnaseAlign,                 
                 histoneAlign = histoneAlign,                
                 histoneNum = histoneNum,                 
                 histoneName = histoneName,                
                 fragL = fragL,                 
                 bowtieInfo = bowtieInfo,                 
                 csemDir = csemDir,                 
                 outfileLoc = outfileLoc,             
                 chrom.ref = chrom.ref))
      
    }else if(!is.null(histoneFile)&is.null(dnaseFile)){
      print("use priorHistone_init and priorHistone_multi functions to process Histone data")
    }else{
      print("Provide Files for deriving Priors or run readAllocate directly without prior information.")
    }
  
}
