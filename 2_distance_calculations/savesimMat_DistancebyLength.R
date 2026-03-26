#### this function creates and saves simMat dissimilarity matrices in csv format for the DistancebyLength.R script 
#### in "~/path/to/Top_Lab/DistancebyLength/[Software/replicate]"
#### for the NN_simMat dissimilarity calculations 

#### parameters:
## data2 <- getdata2(contamination_cutoff,genomesize_cutoff,dir_Similarity_Routput)  ## source("getdata2_DistancebyLength.R") in "~/path/to/Top_Lab/DistancebyLength/
## id_tag <- 'replicate'     ## pick 'replicate' if analysis done per replicate for Software pairs
## id_tag <- 'Software'      ## pick 'Software' if analysis done per software on replicate pairs
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
## complete <- seq(from=0, to=95, by = 5)  ## completeness increments for simMmats
## complete <- c(0,85,90,95)  ## completeness increments for simMmats
## contamination_cutoff <- 10  ## proximeta criteria
## genomesize_cutoff <- c("100Kb"=100e3, "10Mb"=10e6)  ########### genome sizes of bins of 100Kb-10Mb (100e3 to 10e6) based on proximeta criteria,



savesimMat <- function(data2,id_tag,dir_out,complete){
  #### parameters:
  ## data2 <- getdata2(contamination_cutoff,genomesize_cutoff,dir_Similarity_Routput)  ## source("getdata2_DistancebyLength.R") in "~/path/to/Top_Lab/DistancebyLength/
  ## id_tag <- 'replicate'     ## pick 'replicate' if analysis done per replicate for Software pairs
  ## id_tag <- 'Software'      ## pick 'Software' if analysis done per software on replicate pairs
  ## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
  ## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
  ## complete <- seq(from=0, to=95, by = 5) ## completeness increments for simMmats
  ## complete <- c(0,85,90,95)  ## completeness increments for simMmats
  ## contamination_cutoff <- 10  ## proximeta criteria
  ## genomesize_cutoff <- c("100Kb"=100e3, "10Mb"=10e6)  ########### genome sizes of bins of 100Kb-10Mb (100e3 to 10e6) based on proximeta criteria,
  
  
  ##################################################################### id_tag
  if(id_tag=='replicate'){
    pair_tag <- "Software"
    data2$ID <- paste0(as.character(data2$cluster_id),"_",data2[,id_tag],"_GS_",data2$GenomeSize)  ## ex: CL0001_E_15_1 to MetaTOR_122_0_E_15_1
    data2 <- as_tibble(as.data.frame(data2))
    simMat_name <- c("simMat_Proximeta_Bin3C","simMat_Proximeta_MetaTOR", "simMat_Bin3C_MetaTOR")
    simmat_pairs <- str_split_fixed(simMat_name,"_",3)[,-1]
    folders <- unique(data2$replicate)
    objects <- unique(data2$Software)
    print(paste0("id_tag=",id_tag,"; simMat_name=",simMat_name," per ",id_tag))
  }
  if(id_tag=='Software'){
    pair_tag <- "replicate"
    data2$ID <- paste0(as.character(data2$cluster_id),"_",data2[,pair_tag],"_GS_",data2$GenomeSize)  ## ex: CL0001_E_15_1 to CL0001_E_15_2
    data2 <- as_tibble(as.data.frame(data2))
    simMat_name <- c("simMat_E_15_1_E_15_2","simMat_E_15_1_E_15_3A", "simMat_E_15_1_E_15_3B","simMat_E_15_1_E_15_3C",
                     "simMat_E_15_2_E_15_3A","simMat_E_15_2_E_15_3B","simMat_E_15_2_E_15_3C",
                     "simMat_E_15_3A_E_15_3B","simMat_E_15_3A_E_15_3C",
                     "simMat_E_15_3B_E_15_3C")
    tempm <- strsplit(unlist(lapply(unlist(strsplit(simMat_name,"simMat_")), function(x){x[!x ==""]})), split = "_E")
    for (t in 1:length(tempm)){
      tempm[[t]][2] <- paste0("E",tempm[[t]][2])
    }
    simmat_pairs <- matrix(unlist(tempm),nrow=length(tempm),byrow = T)
    folders <- unique(data2$Software)
    objects <- unique(data2$replicate)
    print(paste0("id_tag=",id_tag,"; simMat_name=",simMat_name," per ",id_tag))
  }
  #########################################################################################
  print(paste0("unique(data2$ID):",unique(data2$ID)))
  print("Starting simMats calculations...")  
  for (cc in 1:length(complete)){
    #for (cc in cc){
    print(paste0("Completeness>=",complete[cc]))
    data2_cc = data2 %>% filter(Completeness>=complete[cc])
    
    for (rr in 1:length(folders)){
      #for (rr in r){
      print(paste0("folder=",folders[rr]))
      data2_cc_rr = data2_cc %>% filter(data2_cc[id_tag]==folders[rr])
      contigLengths <- unique(data2_cc_rr[,c("contig","length")])
      nrow(contigLengths)
      length(unique(contigLengths$contig))
      nrow(contigLengths)==length(unique(contigLengths$contig))
      
      ################################  simMat_obj1_obj2 
      for (m in 1:length(simMat_name)){
        simMat_name_m <- simMat_name[m]
        print(paste0("m=",m))
        print(paste0("simMat_name_s=",simMat_name_m))
        
        ###########
        
        fname <- paste0(folders[rr],"_Completeness>=",complete[cc],"_Contamination<=",contamination_cutoff,
                        "_",names(genomesize_cutoff)[1],"-",names(genomesize_cutoff)[2],"_",simMat_name_m)
        if(complete[cc]==0){
          fname <- paste0(folders[rr],"_Completeness>=","ANY","_Contamination<=",contamination_cutoff,
                          "_",names(genomesize_cutoff)[1],"-",names(genomesize_cutoff)[2],"_",simMat_name_m)
        }
        
        path <- paste0(dir_out,folders[rr],"/")
        #setwd(paste0(dir_out,"/",folders[rr]))
        
        if( file.exists(paste0(path,paste0(fname,".csv"))) ){
          #print(paste0("Skipping: ",paste0(fname,".RDS")))
          #print(paste0("Skipping: ",paste0("nnXY_",fname,".RDS")))
          #print(paste0("Skipping: ",paste0("nnYX_",fname,".RDS")))
          print(paste0("Skipping: ",paste0(fname,".csv")))
          #print(paste0("Skipping: ",paste0("nnXY_",fname,".csv")))
          #print(paste0("Skipping: ",paste0("nnYX_",fname,".csv")))
        } else {
          simMat_obj <- c()
          data2_cc_rr_obj1 <- data2_cc_rr %>% filter(data2_cc_rr[pair_tag]==simmat_pairs[m,1])
          data2_cc_rr_obj2 <- data2_cc_rr %>% filter(data2_cc_rr[pair_tag]==simmat_pairs[m,2])
          
          
          ## obj1_obj2
          for(i in unique(data2_cc_rr_obj1$ID)){
            for(j in unique(data2_cc_rr_obj2$ID)){
              s1 <- data2_cc_rr_obj1$contig[data2_cc_rr_obj1$ID == i]
              s2 <- data2_cc_rr_obj2$contig[data2_cc_rr_obj2$ID == j]
              s3 <- intersect(s1,s2)
              s4 <- union(s1,s2)
              num <- sum(contigLengths$length[match(s3,contigLengths$contig)])
              den <- sum(contigLengths$length[match(s4,contigLengths$contig)])
              simMat_obj <- c(simMat_obj, 1 - num/den) #see vegdist: binary version of several distance measures (e.g. canberra, clark, altGower)
              ### distance can more compactly be written as 1 - length(intersect(s1,s2))/length(union(s1,s2))
            }   ############ END for(j in unique(data2_cc_rep_obj2$ID)){
          }  ############ END for(i in unique(data2_cc_rep_obj1$ID)){
          length(unique(data2_cc_rr_obj1$ID))
          length(unique(data2_cc_rr_obj2$ID))
          #dim(simMat_obj)
          
          simMat_obj <- matrix(simMat_obj, nrow = length(unique(data2_cc_rr_obj1$ID)), byrow = T)
          dim(simMat_obj)
          
          rownames(simMat_obj) <- unique(data2_cc_rr_obj1$ID)
          colnames(simMat_obj) <- unique(data2_cc_rr_obj2$ID)
          assign(simMat_name_m,simMat_obj)
          #get(simMat_name_m)
          #saveRDS(object = get(simMat_name_m),file =paste0(fname,".RDS"))
          simMat_df <- as.data.frame(get(simMat_name_m))
          #write.csv(simMat_df,file =paste0(fname,".csv"))
          write.csv(simMat_df,file.path(path,file =paste0(fname,".csv")) )
          #read.csv(file =paste0(fname,".csv"),row.names = 1)
          
          
          #nnYX_simMat_obj<- apply(get(simMat_name_m),1,min) #nnXY is nearest neighbor (highest similarity) to clusters in rep Y from clusters in rep X
          #nnXY_simMat_obj <- apply(get(simMat_name_m),2,min)
          #assign(paste0("nnYX_",simMat_name_m),nnYX_simMat_obj)
          #assign(paste0("nnXY_",simMat_name_m),nnXY_simMat_obj)
          
          
          #saveRDS(object = get(paste0("nnYX_",simMat_name_m)),file =paste0("nnYX_",fname,".RDS"))
          #saveRDS(object = get(paste0("nnXY_",simMat_name_m)),file =paste0("nnXY_",fname,".RDS"))
          #nYX_df <- as.data.frame(get(paste0("nnYX_",simMat_name_m)))
          #nXY_df <- as.data.frame(t(get(paste0("nnXY_",simMat_name_m))))
          #write.csv(nYX_df,file =paste0("nnYX_",fname,".csv"),row.names = T, col.names = F)
          #write.csv(nXY_df,file =paste0("nnXY_",fname,".csv"),row.names = F, col.names = T)
          
        }   ####### END if/else {
      }    ####### END for (m in length(simMat_name)){
    }  ############# END for (rr in length(folders)){
  }   ################ END for (c in length(complete)){
  ######################################################################################################
  #setwd(dir_out)
  return(print(paste0("Analysis complete: simMats saved in [Software/replicate] folders of ",dir_out)))
}



