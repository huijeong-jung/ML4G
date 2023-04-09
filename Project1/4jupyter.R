# Load Libraries

library(dplyr, quietly = T)
library(GenomicRanges, quietly = T)
library(foreach, quietly = T)
library(doParallel, quietly = T)
library(xgboost, quietly = T)
library(ggplot2, quietly = T)

############################ READ IN DATA ####################################

files <- list.files("./ML4G_Project_1_Data/")
bed_files <- files[stringr::str_detect(files,"-bed")]

# Load the beds, join them and convert them to genomic range objects and 
# store them in a list

cell_lines <- c("X1","X2","X3")

beds <- list()
beds2 <- list()

for (cell_line in cell_lines) {
  temp <- list()
  for (bed_file in bed_files) {
    
    tempII <- data.table::fread(paste0(c("./ML4G_Project_1_Data/",bed_file,"/",cell_line,".bed"),collapse = ""))
    colnames(tempII)
    print(tempII)
    if(bed_file == "DNase-bed" ){ 
      tempII$score <- tempII$V7
      tempII$score_scaled <- scale(tempII$V7) + 2
    } else {
      tempII$score <- tempII$V5
      tempII$score_scaled <- scale(tempII$V5) + 2
    }
    
    tempII$type <- bed_file
    
    temp[[bed_file]] <- tempII
    
  }
  
  beds[[cell_line]] <- data.table::rbindlist(temp) 
  beds[[cell_line]]$type <- as.factor(  beds[[cell_line]]$type )
  
  beds2[[cell_line]] <- beds[[cell_line]]
  beds2[[cell_line]]$cell_line <- cell_line
  
  beds[[cell_line]] <-beds[[cell_line]]  %>% GenomicRanges::makeGRangesFromDataFrame(
    keep.extra.columns = T,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3")
  
  
}
