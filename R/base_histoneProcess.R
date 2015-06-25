.histoneProcess=function(histoneFile,fragL,AllocThres,chrList,chrom.ref,capping,outfileLoc="./",outfile,bowtieIndex,csemDir, bowtieDir,vBowtie,mBowtie,pBowtie, save.files){
  if(!file.exists(outfileLoc)){
    system(paste('mkdir ',outfileLoc,sep=''))
  }
  setwd(outfileLoc)
  
  histoneFormat=strsplit(histoneFile,"/")[[1]]
  histoneFormat=histoneFormat[length(histoneFormat)]
  histoneFormat=sub(".*\\.","",histoneFormat)
  # Define filenames and extensions 
  filename=strsplit(histoneFile,"/")[[1]]
  filename=filename[length(filename)]
  filename=strsplit(filename,paste('.',sub(".*\\.","",filename), sep=''))[[1]][1]
  file.sam=paste(filename,'.sam',sep='')
  file.bam=paste(filename,'.bam',sep='')
  file.bed=paste(filename,'.bed',sep='')
  file.bed_AllocThres=paste(filename,'_',AllocThres,'.bed',sep='')
  # Depending on Histone file extension what to do (fastq, bam or bed)
  if(histoneFormat!="bed")
    {
      if(histoneFormat!="bam")
        {
          message( "Info: Aligning reads" )
          FileFormat<-substr(histoneFile, nchar(histoneFile), nchar(histoneFile)+1) 
          if(FileFormat=="a") FileFormat<-"f" 
          ## Generates SAM file
          system(paste(bowtieDir,'/bowtie -',FileFormat,' -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ',
                       bowtieIndex,' ', histoneFile,' ',file.sam, " 2>&1 | tee ", outfileLoc, 
                       "/priorProcessHistoneBowtie_temp.txt", sep=''))  
          system(paste(csemDir,"/run-csem --sam -p ",pBowtie,' ',file.sam,' ' ,fragL,' ',filename,sep='')) 
        }else{
          file.bam = histoneFile
        }
      message( "Info: Allocating multi-reads using CSEM" )
      system(paste(csemDir,"/csem-generate-input --bed ",file.bam,' ',filename,sep='')) 
    }else{
      file.bed = histoneFile
    }
  # Filter out the alignments with posterior probability < AllocThres/1000
  system(paste("awk '$5 >= ",AllocThres," {print $0}' ",file.bed,">", file.bed_AllocThres,sep=""))
  
  message( "Info: Reading the aligned read file and processing it into nucleotide-level files..." )
  .constructNeucleotideLevelCounts(outfileLoc,file.bed_AllocThres,fragL,outfileLoc,chrList,capping)
  countfile=paste('_',file.bed_AllocThres,'_fragL',fragL,'_bin1.txt',sep='')
  
  message( "Info: Partitioning genome..." )
  r=.clusterPositions_3cluster(countfile,chrList,outfileLoc ,chrom.ref,outfile)
  link=vector('list',length(chrList))
  names(link)=chrList
  for(i in chrList){
    link[[i]]=paste(outfileLoc,'/',i,'_',outfile,'_positions_3cluster.txt',sep='')
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
  return(link)
}
