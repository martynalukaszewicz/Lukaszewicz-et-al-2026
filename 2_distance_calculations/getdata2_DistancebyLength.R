#### this function prepares data2 for the DistancebyLength.R script 
#### in "~/path/to/Top_Lab/DistancebyLength"
#### for the simMat pairwise dissimilarity matrix calculations 

#### parameters:
## contamination_cutoff <- 10000
## genomesize_cutoff <- c("100Kb"=100e3, "10Mb"=10e6)
## dir_Similarity_Routput <- "~/path/to/Top_Lab/Similarity_Routput/"
## dir_Similarity_Routput <- "~/path/to/Top_Lab/Similarity_Routput/"



getdata2 <- function(contamination_cutoff,genomesize_cutoff,dir_Similarity_Routput){
  #### parameters:
  ## contamination_cutoff <- 10
  ## genomesize_cutoff <- c("100Kb"=100e3, "10Mb"=10e6)
  ## dir_Similarity_Routput <- "~/path/to/Top_Lab/Similarity_Routput/"
  ## dir_Similarity_Routput <- "~/path/to/Top_Lab/Similarity_Routput/"
  library(readr)
  library(dplyr)
  
  load1 = function(x){read_tsv(x, col_name = TRUE)}
  
  #####################################################################################################
  ###################################################################################### metadata
  ############### read bin3c metadata
  temp = paste0(dir_Similarity_Routput,"bin3c/metadata.tsv")
  bin3c_metadata <- load1(temp) 
  bin3c_metadata <- bin3c_metadata  %>% filter(extent >= as.numeric(genomesize_cutoff[1]) & extent <= as.numeric(genomesize_cutoff[2]))  ## proximeta criteria
  bin3c_metadata <- bin3c_metadata %>% rename("GenomeSize" = "extent")
  bin3c_metadata <- bin3c_metadata  %>% filter(Contamination <= contamination_cutoff)  ## Proximeta and checkm criteria

  ############### read proximeta metadata
  temp = paste0(dir_Similarity_Routput,"proximeta/metadata.tsv")
  proximeta_metadata <- load1(temp) 
  proximeta_metadata <- proximeta_metadata %>% filter(genome_size >= as.numeric(genomesize_cutoff[1]) & genome_size <= as.numeric(genomesize_cutoff[2]))  ## proximeta criteria
  proximeta_metadata <- proximeta_metadata %>% rename("GenomeSize" = "genome_size")
  proximeta_metadata <- proximeta_metadata %>% rename("Contamination" = "marker_gene_overrepresentation")
  proximeta_metadata <- proximeta_metadata %>% rename("Completeness" = "completeness")
  proximeta_metadata <- proximeta_metadata   %>% filter(Contamination <= contamination_cutoff)  ## Proximeta and checkm criteria
  
  ############### read metator metadata
  temp = paste0(dir_Similarity_Routput,"metator/metadata.tsv")
  metator_metadata <- load1(temp) 
  metator_metadata <- metator_metadata  %>% filter(size >= as.numeric(genomesize_cutoff[1]) & size <= as.numeric(genomesize_cutoff[2]))  ## proximeta criteria
  metator_metadata <- metator_metadata %>% rename("GenomeSize" = "size")
  metator_metadata <- metator_metadata  %>% rename("Completeness" = "completness")
  metator_metadata <- metator_metadata  %>% rename("Contamination" = "contamination")
  metator_metadata <- metator_metadata %>% filter(Contamination <= contamination_cutoff)  ## Proximeta and checkm criteria

  
  
  ####################################################################################################
  ###################################################################################### data
  ############### read bin3c data
  temp = paste0(dir_Similarity_Routput,"bin3c/data.tsv")
  bin3c_data <- load1(temp) 
  bin3c_data <- bin3c_data %>% rename("cluster_id" = "cluster")
  bin3c_data$software <- "Bin3C"
  bin3c_data2 <- merge(bin3c_data, bin3c_metadata, by=c("cluster_id","replicate"))
  bin3c_data2 <- bin3c_data2 %>% rename("software" = "software.y")
  bin3c_data2 <- bin3c_data2 %>% rename("Software" = "software.x")
  
  ############### read proximeta data
  temp = paste0(dir_Similarity_Routput,"proximeta/data.tsv")
  proximeta_data <- load1(temp) 
  proximeta_data <- proximeta_data %>% rename("cluster_id" = "cluster")
  proximeta_data$software <- "Proximeta"
  proximeta_data2 <- merge(proximeta_data, proximeta_metadata, by=c("replicate","cluster_id"))
  proximeta_data2 <- proximeta_data2 %>% rename("software" = "software.y")
  proximeta_data2 <- proximeta_data2 %>% rename("Software" = "software.x")

  ############### read metator data
  temp = paste0(dir_Similarity_Routput,"metator/data.tsv")
  metator_data <- load1(temp) 
  metator_data <- metator_data %>% rename("cluster_id" = "cluster")
  metator_data$software <- "MetaTOR"
  metator_data2 <- merge(metator_data, metator_metadata, by=c("replicate","cluster_id"))
  metator_data2 <-merge(metator_data, metator_metadata, by=c("replicate","cluster_id"))
  metator_data2 <- metator_data2 %>% rename("software" = "software.y")
  metator_data2 <- metator_data2 %>% rename("Software" = "software.x")

  #############################################################
  
  
  
  ###################################################################################################################
  ################################################# Create data2:
  ################ combine: replicate, cluster_id, contig, length, software, Software, Completeness, Contamination, GenomeSize
  ############### of proximeta_data2, bin3c_data2, metator_data2 to create data2
  c_names <- c("replicate", "cluster_id", "contig", "length", "software", "Software", "Completeness", "Contamination", "GenomeSize")
  print(paste0("c_names:",c_names))
  #print(paste0("colnames(proximeta_data2):",proximeta_data2))
  print("selecting all_of(c_names)")
  #proximeta_data2_temp <- proximeta_data2 %>% select(all_of(c_names))
  proximeta_data2_temp <- proximeta_data2[,c_names]
  #bin3c_data2_temp <- bin3c_data2 %>% select(all_of(c_names))
  bin3c_data2_temp <- bin3c_data2[,c_names]
  #metator_data2_temp <- metator_data2 %>% select(all_of(c_names))
  metator_data2_temp <- metator_data2[,c_names]
  
  ################################################################## combine data2
  data2 <- rbind(proximeta_data2_temp,bin3c_data2_temp,metator_data2_temp)
  return(data2)
}





