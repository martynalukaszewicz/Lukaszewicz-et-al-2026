############ this function creates and saves NN simMat dissimilarity file called [Software/replicate]_summary.tsv for the DistancebyLength.R script 
#### in "~/path/to/Top_Lab/DistancebyLength/[Software/replicate]"
#### for the plots generation


#### parameters:
## id_tag <- 'replicate'     ## pick 'replicate' if analysis done per replicate for Software pairs
## id_tag <- 'Software'      ## pick 'Software' if analysis done per software on replicate pairs
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"




savesummary <- function(id_tag,dir_out){
  #### parameters:
  ## id_tag <- 'replicate'     ## pick 'replicate' if analysis done per replicate for Software pairs
  ## id_tag <- 'Software'      ## pick 'Software' if analysis done per software on replicate pairs
  ## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
  ## dir_out <- "~/path/to/Top_Lab/DistancebyLength/"
  
  
  ##################################################################### id_tag
  if(id_tag=="replicate"){
    folders <- c("E_15_1",  "E_15_2",  "E_15_3A", "E_15_3B", "E_15_3C")
  }
  if(id_tag=="Software"){
    folders <- c("Proximeta",  "Bin3C",  "MetaTOR")
  }
  
  ######################################################################################
  
  ########################################################### [Software/replicate]_summary.tsv files creation
  for (rr in 1:length(folders)){
    print(paste0("rr=",rr))
    path <- paste0(dir_out,folders[rr],"/")
    print(paste0("folders[rr]=",folders[rr]))
    load2 = function(x){read.csv( x, row.names = 1)}   ################ load all the simMat csv files in that folder
    temp = list.files(path,pattern="*.csv")
    print( paste0("simMat files to summarize in one NN file per:", folders[rr]) )
    ################################# extracting GenomeSize, Completenes cutoffs, Contamination cutoffs, and simMat software pairs 
    ############################### in order to create a single summary table of dissimilarities per folder
    temp_df <- as.data.frame(matrix(unlist(strsplit(temp,split="_Completeness>=")),nrow=length(temp),byrow=T))
    temp_df[,2:3] <- as.data.frame(matrix(unlist(strsplit(temp_df[,2],split="_Contamination<=")),nrow=length(temp),byrow=T))
    temp_df <- cbind(temp_df[,1:2],as.data.frame(matrix(unlist(strsplit(temp_df[,3],split="_")),nrow=length(temp),byrow=T)))
    temp_df <- cbind(temp_df[,-ncol(temp_df)],as.data.frame(matrix(unlist(strsplit(temp_df[,ncol(temp_df)],split=".csv")),nrow=length(temp),byrow=T)))
    
    if(id_tag=="Software"){
      ncol_temp_df <- ncol(temp_df)
      for (nn in 1:ncol_temp_df){
        if (unique(temp_df[,nn]=="simMat")){
          col_idx <- nn
          mid <- 0.5*(ncol_temp_df-col_idx)
        }  ## END if (unique(temp_df[,nn]=="simMat"))
      }  ## END for (nn in 1:ncol_temp_df)
      paste_noNA <- function(x,sep=", ") {
        gsub(pattern=", " ,replacement=sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) ) }
      sep <- "_"
      temp_df$obj1 <- apply( temp_df[,(col_idx+1):(col_idx+mid)], 1, paste_noNA , sep=sep )
      temp_df$obj2 <- apply( temp_df[,(col_idx+mid+1):(ncol_temp_df)], 1, paste_noNA , sep=sep )
      temp_df <- temp_df[,-c( (col_idx+1):(ncol_temp_df)) ]
    }  ## END if(id_tag=="Software")
    
    
    colnames(temp_df) <- c("Folder","Completeness>=","Contamination<=","GenomeSize","simMat","obj1","obj2")
    temp_df$simMat_obj1_obj2 <- paste0(temp_df$simMat,"_",temp_df$obj1,"_",temp_df$obj2)
    temp_df$idx <- 1:nrow(temp_df)
    temp_df["Completeness>="][temp_df["Completeness>="]=="ANY"] <- 0
    temp_df <- as_tibble(temp_df)
    temp_df <-temp_df %>% mutate_at(c("Completeness>=","Contamination<="), as.numeric)
    #temp_df
    
    unique_compl <- unique(temp_df$`Completeness>=`)
    unique_cont <- unique(temp_df$`Contamination<=`)
    unique_gs <- unique(temp_df$GenomeSize)
    
    ######## create an empty summary table per Folder
    unique_simMat_name <- unique(temp_df$simMat_obj1_obj2)
    nnYX_unique_simMat_name <- paste0("nnYX_",unique_simMat_name)
    nnXY_unique_simMat_name <- paste0("nnXY_",unique_simMat_name)
    cluster_id_nnYX_unique_simMat_name <- paste0("cluster_id_",nnYX_unique_simMat_name)
    cluster_id_nnXY_unique_simMat_name <- paste0("cluster_id_",nnXY_unique_simMat_name)
    
    names_r_summary <- c("GenomeSize","Completeness>=","Contamination<=")
    for (nrs in 1:length(unique_simMat_name)){
      names_r_summary <- c(names_r_summary,nnYX_unique_simMat_name[nrs],cluster_id_nnYX_unique_simMat_name[nrs],
                           nnXY_unique_simMat_name[nrs],cluster_id_nnXY_unique_simMat_name[nrs])
    }
    
    ############### list to holder the table information
    r_summary <- sapply(names_r_summary,function(x) NULL)
    nrow_temp_df <- nrow(temp_df)
    
    for (ii in 1:nrow_temp_df){
      print(paste0("ii=",ii))
      print(paste0("temp[ii]=",temp[ii]))
      simMatXY <- load2( paste0(path,temp[ii]) )
      
      nnYX <-apply(simMatXY,1,min) #nnXY is nearest neighbor (highest similarity) of clusters in rep Y to clusters in rep X
      names_nnYX <- names(nnYX)
      nnYX <- nnYX %>% as.vector() 
      
      nnXY <-apply(simMatXY,2,min) 
      names_nnXY <- names(nnXY)
      nnXY <- nnXY %>% as.vector()
      
      len_nnYX <- length(nnYX)
      len_nnXY <- length(nnXY)
      max_len <- max( len_nnYX , len_nnXY )
      min_len <- min( len_nnYX , len_nnXY )
      diff_len <- (max_len-min_len)
      
      #names_r_summary
      #temp_df[s,]
      ######## these always get filled with same values per file simMat file
      r_summary[["GenomeSize"]] <- c(r_summary[["GenomeSize"]], rep(temp_df$GenomeSize[ii],max_len) )
      r_summary[["Completeness>="]] <- c(r_summary[["Completeness>="]], rep(temp_df$`Completeness>=`[ii],max_len) )
      r_summary[["Contamination<="]] <- c(r_summary[["Contamination<="]], rep(temp_df$`Contamination<=`[ii],max_len) )
      
      
      ########################### Dissimilarity
      ####################### nnYX:
      title_temp_nnYX <- paste0("nnYX_",temp_df$simMat_obj1_obj2[ii])
      cluster_id_title_temp_nnYX <- paste0("cluster_id_",title_temp_nnYX)
      r_summary[[cluster_id_title_temp_nnYX]] <- c(r_summary[[cluster_id_title_temp_nnYX]], names_nnYX )
      r_summary[[title_temp_nnYX]] <- c(r_summary[[title_temp_nnYX]], nnYX )
      if (len_nnYX < len_nnXY){
        r_summary[[cluster_id_title_temp_nnYX]] <- c(r_summary[[cluster_id_title_temp_nnYX]], rep(NA,diff_len) )
        r_summary[[title_temp_nnYX]] <- c(r_summary[[title_temp_nnYX]], rep(NA,diff_len) )
      }
      
      ####################### nnXY:
      title_temp_nnXY <- paste0("nnXY_",temp_df$simMat_obj1_obj2[ii])
      cluster_id_title_temp_nnXY <- paste0("cluster_id_",title_temp_nnXY)
      r_summary[[cluster_id_title_temp_nnXY]] <- c(r_summary[[cluster_id_title_temp_nnXY]], names_nnXY )
      r_summary[[title_temp_nnXY]] <- c(r_summary[[title_temp_nnXY]], nnXY )
      if (len_nnXY < len_nnYX){
        r_summary[[cluster_id_title_temp_nnXY]] <- c(r_summary[[cluster_id_title_temp_nnXY]], rep(NA,diff_len) )
        r_summary[[title_temp_nnXY]] <- c(r_summary[[title_temp_nnXY]], rep(NA,diff_len) )
      }
      
      ########################################################## collect indeces of non-filled columns:
      nonfilled_names <- names_r_summary[( (names_r_summary !="GenomeSize") & (names_r_summary !="Completeness>=") & (names_r_summary !="Contamination<=") & 
                                             (names_r_summary !=cluster_id_title_temp_nnYX) & (names_r_summary !=title_temp_nnYX) & 
                                             (names_r_summary !=cluster_id_title_temp_nnXY) & (names_r_summary !=title_temp_nnXY) )]
      
      ############################################ find non-filled variables:
      for(cnames in nonfilled_names ){
        #print(print(paste0("ii=",ii,"; cnames=",cnames, "; temp[ii]=",temp[ii],"; max_len=",max_len)))
        r_summary[[cnames]] <- c(r_summary[[cnames]], rep(NA,max_len) )
      }  ## END for(cnames
      if (ii==nrow_temp_df){
        numeric_cols <- names_r_summary[!str_detect(names_r_summary, "GenomeSize|cluster_id")]
        r_summary <- as.data.frame(do.call(cbind, r_summary))  %>% as_tibble()
        r_summary <- r_summary %>% mutate_at(numeric_cols, as.numeric)
        #setwd(paste0(dir_out,replicates[rr],"/"))
        write_tsv(r_summary,paste0(path,folders[rr],"_summary.tsv"))
      } ## END if (ii==nrow_temp_df)
    }  ## END for (ii
  }  ## for (rr in 1:length(folders))
  ###############################################################################################################################
  return(print(paste0("Analysis complete: [Software/replicate]_summary.tsv saved in [Software/replicate] folders of ",dir_out)))
}



