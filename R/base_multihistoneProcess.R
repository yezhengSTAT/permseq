.multihistoneProcess=function(histoneFile, fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc="./", outfile, bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, save.files){

  #make sure the output folder is created
  if(!file.exists(outfileLoc)){
    system(paste('mkdir ',outfileLoc,sep=''))
  }
  setwd(outfileLoc)

  #find out the input format of histone. BAM/BED input format can save time by skipping the alignment process
  histoneFormat <- strsplit(histoneFile,"/")[[1]]
  histoneFormat <- histoneFormat[length(histoneFormat)]
  histoneFormat <- sub(".*\\.","",histoneFormat)
  
  # Define filenames and extensions
  filename <- strsplit(histoneFile,"/")[[1]]
  filename <- filename[length(filename)]
  filename <- strsplit(filename,paste('.',sub(".*\\.","",filename), sep=''))[[1]][1]
  file.sam <- paste(filename,'.sam',sep='')
  file.bam <- paste(filename,'.bam',sep='')
  file.bed <- paste(filename,'.bed',sep='')
  file.bed_AllocThres <- paste(filename,'_',AllocThres,'.bed',sep='')
  
  # Depending on Histone file extension what to do (fastq, bam or bed)
  if(histoneFormat!="bed")
    {
      if(histoneFormat!="bam")
        {
          #input histone file is in fastq format
          message( "Info: Aligning reads" )
          FileFormat<-substr(histoneFile, nchar(histoneFile), nchar(histoneFile)+1) 
          if(FileFormat=="a") FileFormat<-"f"
          ## Generates SAM file
          system(paste(bowtieDir,'/bowtie -',FileFormat,' -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ',
                       bowtieIndex,' ', histoneFile,' ',file.sam, " 2>&1 | tee ", outfileLoc, 
                       "/priorProcessHistoneBowtie_temp.txt", sep='')) 
          
          ## Get chrList from SAM file
          if(length(chrList)==0)
            { 
              con <- file(file.sam, open = "r") 
              chrList <- c()
              n <- 1
              while(substring(oneLine <- readLines(con, n=1, warn = FALSE), 1, 3) != "@PG"){
                if(grepl("SN", oneLine)){
                  chrList[n] <- gsub('.*SN:(.*)\t.*','\\1', oneLine)
                  n <- n+1
                 } 
              }
            }
          
          system(paste(csemDir,"/run-csem --sam -p ",pBowtie,' ',file.sam,' ',fragL,' ',filename,sep='')) 
        }else{
          file.bam = histoneFile
        }

      message( "Info: Allocating multi-reads using CSEM" )
      system(paste(csemDir,"/csem-generate-input --bed ",file.bam,' ',filename,sep=''))

    }else{
      file.bed = histoneFile
    }

   ## If it wasn't created before or given by the user, creates chrList from the bed file
  if(length(chrList)==0){
    system(paste("awk '{print $1}' ",file.bed," | sort | uniq > ",outfileLoc,"/chrList_temp",sep=""))
    chrList=read.table(paste(outfileLoc,"/chrList_temp",sep=""))$V1
    system(paste("rm -rf ",outfileLoc,"/chrList_temp", sep=""))
  }

  # Filter out the alignments with posterior probability < AllocThres/1000
  system(paste("awk '$5 >= ",AllocThres," {print $0}' ",file.bed,">", file.bed_AllocThres,sep=""))
  
  message( "Info: Reading the aligned read file and processing it into nucleotide-level files..." )
  .constructNeucleotideLevelCounts(outfileLoc,file.bed_AllocThres,fragL,outfileLoc,chrList,capping)
  countfile=paste('_',file.bed_AllocThres,'_fragL',fragL,'_bin1.txt',sep='')
  
  message( "Info: Partitioning genome..." )
  r=.clusterPositions_histoneonly(countfile,chrList,outfileLoc ,chrom.ref,outfile)
  dnaseKnots=r[[1]]
  dnaseThres=r[[2]]

  link2=link=vector('list',length(chrList))
  names(link2)=names(link)=chrList
  for(i in chrList){
    link[[i]]=paste(outfileLoc,'/',i,'_',outfile,'_positions_cluster.txt',sep='')
    link2[[i]]=paste(outfileLoc,'/',i,'_',outfile,'_positions_3cluster.txt',sep='')
  }

  # Remove files if save.files=FALSE
  if(save.files=='FALSE' | save.files=='F')
  {
    system(paste("rm -rf ",file.bam,sep = ""))
    system(paste("rm -rf ",file.sam,sep = ""))
    system(paste("rm -rf *sorted.bam*",sep = ""))
    system(paste("rm -rf ",'*_',AllocThres,'.bed*',sep = ""))
    system(paste("rm -rf *_output_Dnasequantile.txt*",sep = ""))
  }

  priorObject=list(dnaseKnots=dnaseKnots,dnaseThres=dnaseThres,posLoc_bychr=link,posLoc3_bychr=link2)
  return(priorObject)
}
