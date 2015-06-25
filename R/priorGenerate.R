priorGenerate=function(object,chipFile,maxHistone=6,outfileLoc="./"){

  if(length(chipFile)>1){
    message( "More than one ChIP-seq data detected: Replicates must be from the same TF (experiment). Otherwise, we recommend run the different ChIP-seq data separately.")
  }
  if(!file.exists(outfileLoc)){
    system(paste('mkdir ',outfileLoc,sep=''))
  }

  # Obtain information from "object"
  bowtieIndex=object['bowtieInfo']$bowtieIndex
  bowtieDir=object['bowtieInfo']$bowtieDir
  vBowtie=object['bowtieInfo']$vBowtie
  mBowtie=object['bowtieInfo']$mBowtie
  pBowtie=object['bowtieInfo']$pBowtie
  csemDir=object['csemDir']
  chrlist=names(object['posLoc_bychr'])
  dnaseKnots=object['dnaseKnots']
  dnaseThres=object['dnaseThres']
  object['chipNum'] = length(chipFile)

  setwd(outfileLoc)
  outfile=NULL
  # Name of ChIP file(s)
  for(i in 1:length(chipFile)){
    t1=strsplit(chipFile[i],'/',)[[1]]
    t1=t1[length(t1)]
    t=strsplit(t1,paste('.',sub(".*\\.","",t1), sep=''))[[1]][1]
    outfile=c(outfile,t)
  }
  
  if(length(object['chipSAM'])==0){ # If ChIP data hasn't been processed before
    message( "Info: Aligning reads using BOWTIE..." )
    chipAlign = list()
    for(i in 1:length(chipFile)){
      FileFormat<-substr(chipFile[i], nchar(chipFile[i]), nchar(chipFile[i])+1)
      if(FileFormat=="a"){
        FileFormat<-"f"
      }
      system(paste(bowtieDir,'/bowtie -',FileFormat,' -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ', bowtieIndex, ' ',chipFile[i], ' ',outfile[i],'.sam', " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", outfile[i], "_Bowtie_temp.txt",sep=''))
      
      chipAlign[[outfile[i]]] = .summary(outfileLoc, paste("priorProcessChIP_", outfile[i], "_Bowtie_temp.txt", sep = ""))
    }

    object['chipAlign'] = chipAlign
    object['chipSAM']=paste(outfileLoc,'/',outfile,'.sam',sep='')
  }
  outfile_chipmean=paste(outfile,'_chipmean_temp',sep='')
  
  message( "Info: Calculating averaged ChIP counts..." )
  object=.chipMeanCounts_multi(object, outfileLoc, outfile_chipmean)
  
  message( "Info: Variable selection..." )
  num_require=2+length(dnaseKnots)+maxHistone*2
  if(object['dataNum']>1){ 
    setwd(outfileLoc)
    index=c(rep(1,4),NA,rep(2:object['dataNum'],each=2))
    # Group Lasso
    result=mclapply(outfile_chipmean,function(x){.run_grplasso(x,dnaseKnots,dnaseThres,chrlist,index,num_require)},mc.preschedule = FALSE)
    choose=data_all=vector('list',length(outfile_chipmean))

    for(i in 1:length(outfile_chipmean)){
      # For each ChIP data - choose: selected histone(s) data, data_all_combined: aggregated data
      choose[[i]]=result[[i]][['choose']]
      data_all[[i]]=result[[i]][['data_all_combined']]
    }
    
    inter=filter=vector('list',min(unlist(lapply(choose,length))))
    k=min(c(object['dataNum']-1,maxHistone))
    inter=choose[[1]][1:k]
    if(length(object['chipSAM'])>1){
      for(j in 2:(length(object['chipSAM'])))
        inter=unique(intersect(inter,choose[[j]][1:k]))
    }
    inter=inter[order(as.numeric(inter),decreasing=F)]
    if(length(inter)!=0){
      filter=c(1:(object['dataNum']-1))[-inter]
      }else{
        filter=c(1:(object['dataNum']-1))
        inter=NULL
      }
    if(length(filter)==0) filter=NULL
    }else{
      data_all=vector('list',length(outfile_chipmean))
      for(i in 1:length(outfile_chipmean)){
        file=outfile_chipmean[i]
        for(j in 1:length(chrlist)){
          chr.i=chrlist[j]
          data_all[[i]]=rbind(data_all[[i]],t(.fast.read.table(paste(outfileLoc,'/',chr.i,'_',file,sep=''),3)))
        }
      data_all[[i]]=.aggregate_data(data_all[[i]])
    }
    filter=NULL
    inter=NULL
  }

  message( "Info: Model fitting..." )
  y=vector('list',length(chipFile))

  for(i in 1:length(chipFile)){
    y[[i]]=.fitRlm(data_all[[i]],filter,object)
  }

  message( "Info: Generating prior files..." )
  setwd(outfileLoc)
  CMD='cat '
  for(i in 1:length(object['posLoc_bychr'])){
    CMD=paste(CMD,object['posLoc_bychr'][[i]],sep=' ')
  }
  CMD=paste(CMD,'>cluster_temp.txt',sep='')
  system(CMD)

  positionFileLoc_tfs='cluster_temp.txt'
  tt=0:(length(object['dnaseThres'])-1)
  tt=tt[order(as.numeric(tt),decreasing=F)]
  if(object['dataNum']>1){
    data_id=t(apply(as.matrix(data_all[[1]][,1]),1,function(x){strsplit(as.character(x),'_')[[1]]}))
    tt=data_id[,1]
    for(i in inter){
      tt=paste(tt,data_id[,1+i],sep='_')
    }
    tt=unique(tt)
    tt=matrix(as.character(tt),nrow=1)
    write.table(tt,file='dnase_pattern_id_temp.txt',quote=F,col.names=F,row.names=F,sep='\t')
    script='select_multi_pattern.pl'
    Fn.Path <- system.file( file.path("Perl",script), package="permseq")
    CMD<-paste('perl  ',Fn.Path,'dnase_pattern_id_temp.txt','cluster_temp2.txt',positionFileLoc_tfs,sep=' ')
    for(i in c(0,inter)){
      CMD=paste(CMD,i,sep=' ')
    }
    system(CMD)
    positionFileLoc_tfs='cluster_temp2.txt'
  }
  for(i in 1:length(chipFile)){
    parameter_a=c(length(tt),as.numeric(y[[i]][match(tt,y[[i]][,1]),2])+1) 
    write.table(matrix(parameter_a,nrow=1),file='rep_dnase_parameter_temp.txt',col.names=F,row.names=F,quote=F)
    system(paste('cat rep_dnase_parameter_temp.txt ',positionFileLoc_tfs,'>prior_',outfile[i],'.txt',sep=''))
  }
  system('rm -rf *temp*')
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
  if(is.null(inter)){
    object['prior'] <- paste(outfileLoc,'/','prior_',outfile,'.txt',sep='')
    object['outfileLoc'] <- outfileLoc
    object['chipName'] <- outfile
    return(object)
  }else{
    object['prior'] <- paste(outfileLoc,'/','prior_',outfile,'.txt',sep='')
    object['outfileLoc'] <- outfileLoc
    object['histoneGrpL'] <- object['histoneName'][inter]
    object['chipName'] <- outfile
    return(object)
  }
}
