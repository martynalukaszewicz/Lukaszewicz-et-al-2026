############ this function creates and saves removed NA's form NN simMat dissimilarity file called NN_[Software/replicate]_summary.tsv for the DistancebyLength.R script 
#### in "~/path/to/Top_Lab/DistancebyLength/[Software/replicate]"
#### for the plots generation


#### parameters:
## id_tag <- 'replicate'     ## pick 'replicate' if analysis done per replicate for Software pairs
## id_tag <- 'Software'      ## pick 'Software' if analysis done per software on replicate pairs
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"



saveNNsummary <- function(id_tag,dir_out){
  ##################################################################### id_tag
  if(id_tag=="replicate"){
    folders <- c("E_15_1",  "E_15_2",  "E_15_3A", "E_15_3B", "E_15_3C")
  }
  if(id_tag=="Software"){
    folders <- c("Proximeta",  "Bin3C",  "MetaTOR")
  }
  ######################################################################################
  NN_names_r_summary <- c("GenomeSize","Completeness>=","Contamination<=","NN Value","cluster_id_NN","NN_name","NN Pairs")
  ######################################################################## NN_[Softaware/replicate]_summary.tsv
  for (rr in 1:length(folders)){
    NN_r_summary <- sapply(NN_names_r_summary,function(x) NULL)
    
    ############## set working directory to replicate folder
    path <- paste0(dir_out,folders[rr],"/")
    r_summary <- read_tsv(paste0(path,folders[rr],"_summary.tsv"))
    names_r_summary <- names(r_summary)
    unique_gs <- unique(r_summary$GenomeSize)
    unique_compl <- unique(r_summary$`Completeness>=`)
    unique_cont <- unique(r_summary$`Contamination<=`)
    
    print(paste0("rr=",rr))
    print(paste0("folder=",folders[rr]))
    
    for (gs in unique_gs){
      print(paste0("gs=",gs))
      for (compl in unique_compl){
        print(paste0("compl=",compl))
        for (cont in unique_cont){
          print(paste0("cont=",cont))
          r_summary_gs_compl_cont <- r_summary %>% filter( (GenomeSize == gs) & (`Completeness>=` == compl) & (`Contamination<=` == cont) )
          r_summary_gs_compl_cont
          cluster_id_cols <- names_r_summary %>% str_subset(pattern = "cluster_id")
          nn_cols <- gsub("cluster_id_","",cluster_id_cols)
          for (cc in 1:length(nn_cols)){  ## remove NA's from the columns
            temp_cluster_id_NN <- r_summary_gs_compl_cont[cluster_id_cols[cc]]
            length_NN <- length(temp_cluster_id_NN[!is.na(temp_cluster_id_NN)])
            
            NN_r_summary$cluster_id_NN <- c( NN_r_summary$cluster_id_NN, 
                                             (temp_cluster_id_NN[!is.na(temp_cluster_id_NN)]) )  ## fill cluster_id_NN
            
            temp_NN_Value <- r_summary_gs_compl_cont[nn_cols[cc]]
            
            NN_r_summary$`NN Value` <- c(NN_r_summary$`NN Value`, 
                                         (temp_NN_Value[!is.na(temp_NN_Value)]) ) ## fill NN Value
            
            NN_r_summary$GenomeSize <- c( NN_r_summary$GenomeSize, rep(gs,length_NN) )   ## fill GenomeSize
            NN_r_summary$`Completeness>=` <- c( NN_r_summary$`Completeness>=`, rep(compl,length_NN) )   ## fill Completeness
            NN_r_summary$`Contamination<=` <- c( NN_r_summary$`Contamination<=`, rep(cont,length_NN) ) ## fill Contamination
            
            NN_r_summary$NN_name <- c( NN_r_summary$NN_name, rep(nn_cols[cc],length_NN) )  ## fill 
            temp_NN <- unlist(strsplit(unlist(strsplit(nn_cols[cc],split="_simMat_")),"_"))
            ##################################################################### id_tag
            if(id_tag=="Software"){
              temp_NN <- c(temp_NN[1],gsub(pattern=" ",replacement="_",paste(temp_NN[2:4], collapse = " ")),
                           gsub(pattern=" ",replacement="_",paste(temp_NN[5:7], collapse = " ")) )
            }
            ###################################################################
            if("nnYX" %in% temp_NN){
              NN_title <- paste0(temp_NN[3]," to ", temp_NN[2])
            } else {
              NN_title <- paste0(temp_NN[2]," to ", temp_NN[3])
            }
            NN_r_summary$`NN Pairs` <- c(NN_r_summary$`NN Pairs`, rep(NN_title,length_NN) )
          }  ## END for (cc in 1:length(nn_cols))
        }  ## END or (cont in unique_cont)
      }  ## END for (compl in unique_compl)
    }  ## END for (gs in unique_gs)
    NN_r_summary <- as.data.frame(do.call(cbind, NN_r_summary))  %>% as_tibble()
    NN_r_summary <- NN_r_summary %>% mutate_at(c("Completeness>=","Contamination<=","NN Value"), as.numeric)
    write_tsv(NN_r_summary,paste0(path,"NN_",folders[rr],"_summary.tsv"))
  } ## END for (rr in 1:length(folders))  
  return(print(paste0("Analysis complete: NN_[Software/replicate]_summary.tsv saved in [Software/replicate] folders of ",dir_out)))
}



