# Load Libs
library(dplyr)
library(GenomicRanges)
library(foreach)
library(doParallel)
library(xgboost)
library(caret)
library(ggplot2)

############################ READ IN DATA ####################################
genesTrain <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X1_train_info.tsv")
genesVal <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X1_val_info.tsv")
genesTest <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X3_test_info.tsv")

genes <- rbind(genesTrain,genesVal) %>% 
  rbind(genesTest) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    keep.extra.columns = T,
    seqnames.field = "chr",
    start.field = "TSS_start",
    end.field = "TSS_end",
    strand.field = "strand")


# define gene set indices store in list for for loop

train_and_val_inds <- which(genes$gene_name %in% c(genesTrain$gene_name,genesVal$gene_name))
test_inds <- which(genes$gene_name %in% genesTest$gene_name)

cell_line_gene_inds <- list(train_and_val_inds,
                         train_and_val_inds,
                         test_inds)


rF <- 5000

# Use the promoters expansion as the ranges are the TSS

shortIntevals <- promoters(genes,
                          upstream  = rF,
                          downstream= rF)

longIntervals <- promoters(genes,
                          upstream  = 500000 ,
                          downstream= 500000)

geneOverlaps <- findOverlaps(longIntervals, shortIntevals) %>%
                as.data.frame()

geneOverlaps <- split(geneOverlaps, f=geneOverlaps$queryHits)


# Register the parallel backend

num_cores <- 3
registerDoParallel(cores = num_cores)


########################## CREATE FEATURES ##################################

res <- foreach(j = 1:3, .packages = c("dplyr", "GenomicRanges")) %dopar% {
  features  <- list() # Initialize features for cell line j
  
  shortOverlaps <- findOverlaps(shortIntevals, beds[[cell_lines[j]]]) %>%
                   as.data.frame() %>%
                   mutate(queryHits = factor(queryHits,levels = 1:length(genes))) %>%
                   group_by(queryHits,.drop = FALSE) %>%
                   group_split()
  
  LongOverlaps <- findOverlaps(longIntervals, beds[[cell_lines[j]]]) %>%
                  as.data.frame() %>%
                  mutate(queryHits = factor(queryHits,levels = 1:length(genes))) %>%
                  group_by(queryHits,.drop = FALSE) %>%
                  group_split()
            
  
  for (i in cell_line_gene_inds[[j]]) {
    
    feature_bed <- numeric(178) # Number of features
    
    shortInteval <- shortIntevals[i]
          
    longInterval <- longIntervals[i]
    
    geneNeigbourhoodData <- shortIntevals[geneOverlaps[[i]]$subjectHits]
    
    shortRangeData <- beds[[cell_lines[j]]][shortOverlaps[[i]]$subjectHits] 
    
    longRangeData <- beds[[cell_lines[j]]][LongOverlaps[[i]]$subjectHits]
    
    peaks_close_to_a_tssInds <- findOverlaps(geneNeigbourhoodData,longRangeData)
    peaks_not_close_to_a_tssData <- longRangeData[-(peaks_close_to_a_tssInds@to)]
    
    offset <- 0
    for (bed_file in bed_files) {
      
      peaks <- shortRangeData[shortRangeData$type ==  bed_file]
      peaks_distal <- longRangeData[longRangeData$type ==  bed_file]
      peaks_lonely_distal <- peaks_not_close_to_a_tssData[peaks_not_close_to_a_tssData$type ==  bed_file]
      
      feature_bed[1+offset] <- length(peaks) #number_of_peaks
      feature_bed[2+offset] <- length(peaks_distal) - length(peaks) #number_of_distal_peaks
      feature_bed[3+offset] <- length(peaks_lonely_distal) # number_of_distal_lonely_peaks

      
      if (feature_bed[1+offset] > 0) {


        peak_centers <- (start(peaks) + end(peaks)) / 2
        
        center_dis <- peak_centers - ((start(shortInteval) + end(shortInteval))) / 2
        feature_bed[4+offset] <- center_dis[which.min(abs(center_dis))] #distance_to_closest_peak_signed
        feature_bed[5+offset] <- center_dis[which.max(abs(center_dis))] #distance_to_furthest_peak_signed
        feature_bed[6+offset] <- min(abs(center_dis))                 # distance_to_closest_peak_unsigned
        feature_bed[7+offset] <- sum(peaks$score) #sum_of_scores
        feature_bed[8+offset] <- peaks[which.min(abs(center_dis)), ]$score #score_of_closet_peak
        feature_bed[9+offset] <- peaks[which.max(abs(center_dis)), ]$score #score_of_furthest_peak
        feature_bed[10+offset] <- max(peaks$score) #max_peak
        feature_bed[11+offset] <- mean(peaks$score) #mean_peak
        
        # Binned Peak score sum vectors
        
        intervals <- seq(start(shortInteval),end(shortInteval) + 50,by=1000)
        intervals[1]  <- intervals[1]  - 100000
        intervals[11] <- intervals[11] + 100000
        bincodes <- .bincode(peak_centers,intervals)
        for (r in 1:length(peaks)) {
          feature_bed[11+offset+bincodes[r]] <- feature_bed[11+offset+bincodes[r]] + peaks[r]$score
        }
        
        
      } else{
        
        feature_bed[4+offset] <- 1000000 #distance_to_closest_peak_signed
        feature_bed[5+offset] <- 1000000 #distance_to_furthest_peak_signed
        feature_bed[6+offset] <- 1000000 #distance_to_closest_peak_unsigned
        feature_bed[7+offset] <- 0 #sum_of_scores
        feature_bed[8+offset] <- 0 #score_of_closet_peak
        feature_bed[9+offset] <- 0 #score_of_furthest_peak
        feature_bed[10+offset] <- 0 #max_peak
        feature_bed[11+offset] <- 0#mean_peak
        
      }
      if (feature_bed[2+offset] > 0) {
        feature_bed[22+offset] <- sum(peaks_distal$score) - feature_bed[7+offset] # sum_of_scores_distal
        feature_bed[23+offset] <- feature_bed[22+offset]/feature_bed[2+offset] # mean_of_score_distal
        # Add distance to closest distal peak ??? 
      } else {
        feature_bed[22+offset] <- 0 # sum_of_scores_distal
        feature_bed[23+offset] <- 0 # mean_of_score_distal
      }
      if(feature_bed[3+offset]> 0){
        feature_bed[24+offset] <- max(peaks_lonely_distal$score) # lonely_distal_max
        feature_bed[25+offset] <- mean(peaks_lonely_distal$score) # lonely_distal_mean
      } else {
        feature_bed[24+offset] <- 0 # lonely_distal_max
        feature_bed[25+offset] <- 0 # lonely_distal_mean
      }
      offset <- offset + 25
    }


    feature_bed[176] <- (genes[i]$gene_end - genes[i]$gene_start) # gene_length
    feature_bed[177] <- length(geneNeigbourhoodData) - 1 # number_of_genes_in_large_neighbourhood
    feature_bed[178] <-  as.numeric(strand(genes[i]) == "+") # Strand
    
    
    features[[genes[i]$gene_name]] <- feature_bed
    
  }
  features
}
stopImplicitCluster()

# load response vars

train_vals_x1 <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X1_train_y.tsv")
train_vals_x2 <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X2_train_y.tsv")

validation_vals_x1 <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X1_val_y.tsv")
validation_vals_x2 <- data.table::fread("ML4G_Project_1_Data/CAGE-train/X2_val_y.tsv")

features <- c("number_of_peaks",
              "number_of_distal_peaks",
              "number_of_distal_lonely_peaks",
              "distance_to_closest_peak_signed",
              "distance_to_furthest_peak_signed",
              "distance_to_closest_peak_unsigned",
              "sum_of_scores",
              "score_of_closet_peak",
              "score_of_furthest_peak",
              "max_peak",
              "mean_peak",
              paste0("bin_",1:10),
              "sum_of_scores_distal",
              "mean_of_score_distal",
              "lonely_distal_max",
              "lonely_distal_mean"
)


feature_cols <- c()
for (bed_file in bed_files) {
  feature_cols <- c(feature_cols,paste0(features,"_",strsplit(bed_file,"-bed")[[1]]))
}
feature_cols <- c(feature_cols,"gene_length","number_of_genes_in_large_neighbourhood","strand")

# GET DATA INTO RIGHT ORDER

examine <- data.frame(res[[1]]) 
rownames(examine) <- feature_cols


x1_feat <- data.frame(t(data.frame(res[[1]])))
colnames(x1_feat) <- feature_cols

rownames(x1_feat) <-  c(train_vals_x1$gene_name,validation_vals_x1$gene_name)
x1_feat$y_val <- c(train_vals_x1$gex,validation_vals_x1$gex)


x2_feat <- data.frame(t(data.frame(res[[2]])))
colnames(x2_feat) <- feature_cols

rownames(x2_feat) <-  c(train_vals_x2$gene_name,validation_vals_x2$gene_name)
x2_feat$y_val <- c(train_vals_x2$gex,validation_vals_x2$gex)

x3_feat <- data.frame(t(data.frame(res[[3]])))
colnames(x3_feat) <- feature_cols
rownames(x3_feat) <-  c(genesTest$gene_name)


train_complete <-data.frame(rbind(x1_feat[train_vals_x1$gene_name,],
                                  x2_feat[train_vals_x2$gene_name,]))

val_complete <-data.frame(rbind(x1_feat[rownames(x1_feat) %in% validation_vals_x1$gene_name  ,],
                                x2_feat[rownames(x2_feat) %in% validation_vals_x2$gene_name,]))

test_complete <- x3_feat

errything <- rbind(list(train_complete,val_complete,test_complete))

#make this example reproducible
set.seed(0)

# Test using only certain features
train_complete<- train_complete %>%
  select(-contains("H3K9me3"))

val_complete<- val_complete %>%
  select(-contains("H3K9me3"))

train_complete<- train_complete %>%
  select(-contains("H3K9me3"))


#define final training and testing sets
 

xgb_train = xgb.DMatrix(data =  data.matrix(train_complete[,-ncol(train_complete)]),
                        label = train_complete[,ncol(train_complete)])

xgb_test = xgb.DMatrix(data = data.matrix(val_complete[,-ncol(val_complete)]),
                       label = val_complete[,ncol(val_complete)])


xgb_real_test = xgb.DMatrix(data = data.matrix(test_complete))


watchlist = list(train=xgb_train, test=xgb_test)


evalerror <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  err <- cor(labels , preds,method = "spearman")
  return(list(metric = "spearman", value = err))
}

#fit XGBoost model and display training and testing data at each round

param <- list(objective = "count:poisson",
              eval_metric = evalerror,
              subsample = 0.7, # does best between 0.6 -0.8
              lambda = 1, #raising helps but not in cuncjuntion with alpha
              alpha=1000, # helps test loss a lot
              colsample_bytree = 0.5 # helps a bit
              )

model = xgb.train(data = xgb_train,
                  max.depth = 7,
                  watchlist=watchlist,
                  nrounds = 500,
                  params = list(objective = "count:poisson",
                            eval_metric = evalerror,
                           subsample = 0.7,
                           colsample_bytree= 0.4,
                           lambda = 1000,
                           alpha = 0))



############################ PARAM SEARCH ####################################


test_pearsons <- list()

# best test pearsons (0.776791) seems to be: max:7 subsample:0.8
# best was  5 500 0.4 1 with 0.7811712
# new best 6 750 0.4 1000 with 0.7817005
for (max_depth in 4:7) {
  for (alpha in c(1,250,500,750,1000)) {
    for (colsample in c(0.3,0.4,0.5)) {
      for (lambda in c(1,500,1000)) {
  
      
      param <- list(objective = "count:poisson",
                    eval_metric = evalerror,
                    subsample = 0.7,
                    colsample_bytree= colsample,
                    lambda = lambda,
                    alpha = alpha)
      
        model = xgb.train(data = xgb_train,
                          max.depth = max_depth,
                          watchlist=watchlist,
                          nrounds = 500,
                          params = param,
                          verbose = 0)
        test_pearsons[[paste(max_depth,
                             alpha,
                             colsample,
                             lambda)]] <- model$evaluation_log$test_spearman
      }
    }
  }
}


losses_df <- data.frame(
  iteration = 1:length(test_pearsons[[1]]),
  stack(sapply(test_pearsons, `length<-`, max(lengths(test_pearsons))))
)

# Rename the columns
colnames(losses_df) <- c("iteration", "iter", "hyperparameter","loss")

# Plot the losses using ggplot2
ggplot(data = losses_df, aes(x = iteration, y = loss, color = as.factor(hyperparameter))) +
  geom_line() +
  labs(title = "Losses for Different Hyperparameters",
       x = "Iteration",
       y = "Pearsons on Validation Set",
       color = "Max Depth") +
  theme_minimal()



for (i in seq(150,400,25)) {

  final = xgboost(data = xgb_train,
                  max.depth = 6,
                  nrounds = i,
                  params = list(objective = "count:poisson",
                                subsample = 0.7,
                                 colsample_bytree= 0.4,
                                 lambda = 1000,
                                 alpha = 750),
                  verbose = 0)
  
  pred_y = predict(final, xgb_test)
  print(i)
  print(cor(val_complete[,ncol(val_complete)] , pred_y,method = "spearman"))
}
  


final = xgboost(data = xgb_train,
                max.depth = 6,
                nrounds = 250,
                params = list(objective = "count:poisson",
                              subsample = 0.7,
                              colsample_bytree= 0.4,
                              lambda = 1000,
                              alpha = 750))

feature_importance <- xgboost::xgb.importance(model = final)
xgboost::xgb.plot.importance(feature_importance)
xgb.plot.shap(data.matrix(train_complete[,-ncol(train_complete)]), model = final)
xgb.plot.shap(data.matrix(train_complete[,-ncol(train_complete)]), model = final,top_n = 10)
xgb.plot.shap.summary(data.matrix(train_complete[,-ncol(train_complete)]), model = final)


pred_y = predict(final, xgb_real_test)

preds <- data.frame("gene_name" = genesTest$gene_name,
                    "gex_predicted" = pred_y,
                    row.names = 0:(length(pred_y)-1))

write.csv(preds,file = "gex_predicted.csv",row.names = T,quote = FALSE)

lm1 <- train(y~., data = train, method = "lm")
lm1_Preds <- predict(lm1,test)
cor(test_y, lm1_Preds,method = "spearman")


rf1 <- train(y~., data = train, method = "rf")
rf1_Preds <- predict(rf1,test)
cor(test_y, rf1_Preds,method = "spearman")

