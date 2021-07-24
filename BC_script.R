#-----------------------------------------------------------------------------------
# Prediction of molecular subtypes of Breast Cancer based on TCGA proteomic data #
#-----------------------------------------------------------------------------------

## Data modification ##

### loading libraries
library(stabs)
library(factoextra)
library(NbClust)
library(FunCluster)
library(ggfortify)
library(glmnet)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(randomForest)
library(data.table)
library(mlr)
library(h2o)
library(caret)
library(plsVarSel)
library(pROC)
library(jtools)
library(dplyr)
library(cluster)
library("flexclust")
library(corrplot)
library(tsne)
library(clusterCrit)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(randomForest)
library(Rtsne)
library(neuralnet)

###loading the data
clinical <- read.csv("BC_data/clinical_data_breast_cancer.csv") 
genes <- read.csv("BC_data/PAM50_proteins.csv")
proteomes <- read.csv("BC_data/77_cancer_proteomes_CPTAC_itraq.csv")

###proteomes data modification
n <- proteomes$RefSeq_accession_number
proteomes <- as.data.frame(t(proteomes[,4:86]))
colnames(proteomes) <- n #add refseq
proteomes <- cbind(rownames(proteomes), data.frame(proteomes, row.names=NULL)) #separate the sample ID
colnames(proteomes)[1] <- "Complete.TCGA.ID" #the same name of column in clinical data

###sample ID
get_sample_id <- function(proteome_id) {
  x = substr(proteome_id, 4, 7)
  y = substr(proteome_id, 0, 2)
  paste("TCGA",y,x,sep="-")
}
proteomes$Complete.TCGA.ID <- sapply(proteomes$Complete.TCGA.ID, get_sample_id)
proteomes_id <- proteomes

###missing data
NA_counts <- colSums(is.na(proteomes))/nrow(proteomes)
pdf("1. Data modification/BC_missingData.pdf") #print plot to pdf
print(plot(sort(NA_counts, decreasing = TRUE), col = 'blue', type = 'h', xlab = "index of proteome", 
     ylab="proportion of missing data", main = "Missing data for proteome data"))
dev.off()
length(NA_counts[NA_counts>0.25]) #2251 missing data

#remove variables with >25% missing data
proteomes <- proteomes[ , colSums(is.na(proteomes))  / nrow(proteomes) < 0.25]
for (i in which(sapply(proteomes, is.numeric))) { #loop to impute means for remaining missing data
  proteomes[is.na(proteomes[, i]), i] <- mean(proteomes[, i],  na.rm = TRUE)
}
data <-  inner_join(clinical, proteomes, by = "Complete.TCGA.ID") #inner join to create full data
colnames(data)[3] <- "diag_age" #col name change for 3 column

###barplot
ggplot(data, aes(x =data$PAM50.mRNA, fill=data$PAM50.mRNA)) + geom_bar() + ggtitle("Breast Cancer Subtype") 
ggsave('1. Data modification/BC_subtypes_barplot.pdf')

#-----------------------------------------------------------------------------------

## PCA and T-SNE Analysis ##

###PCA
data.pca <- prcomp(data[,32:length(data)], center = TRUE, scale. = TRUE)
summary(data.pca)
autoplot(data.pca, data = data, colour = 'PAM50.mRNA') + ggtitle("PCA Analysis for all BC subtypes")
ggsave('2. PCA and T-SNE Analysis/BC_subtypes_PCA.pdf')

###T-SNE
data.tsne <- data
data.tsne_unique <- unique(data.tsne)
data.tsne <- as.matrix(data.tsne_unique[,32:length(data.tsne)])
set.seed(1)
data.tsne_out <- Rtsne(data.tsne, pca=FALSE, perplexity=1,dims=2, theta=0.5) # Run TSNE
summary(data.tsne_out$Y)
tsne_plot <- data.frame(x = data.tsne_out$Y[,1], y = data.tsne_out$Y[,2], col = data.tsne_unique$PAM50.mRNA) #TSNE plot
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + ggtitle("T-SNE Analysis for all BC subtypes")
ggsave('2. PCA and T-SNE Analysis/BC_subtypes_T-SNE.pdf')

#-----------------------------------------------------------------------------------

## Feature selection ##

###function to repeat lasso regression and return the selected model variables
lasso_sub=function(k=1, x_data, y_data) {
  set.seed(k)
  s=sample(nrow(data), size = 0.8*nrow(data))
  x_sub=x_data[s,]
  y_sub=y_data[s]
  model.sub=cv.glmnet(x=x_sub, y=y_sub, alpha=1, family='multinomial')
  coef.sub=coef(model.sub, s='lambda.1se')[-1]
  return(coef.sub)
}

options(warn = -1) #turn off warnings

###run model 100 times
niter=100
lasso.stab=sapply(1:niter, FUN = lasso_sub, x_data = as.matrix(data[,31:ncol(data)]), 
                  y_data=as.matrix(data[,21]))

###create a matrix of all predictor variables
lasso_matrix <- matrix(nrow=length(lasso.stab[[1]]),ncol=length(lasso.stab))
rownames(lasso_matrix) <- rownames(lasso.stab[[1]])
for (i in 1:300){ #loop through to put list contents into matrix
  lasso.data.frame <- as.matrix(lasso.stab[[i]])
  lasso_matrix[,i] <- lasso.data.frame
}
lasso_matrix <- ifelse(lasso_matrix != 0, 1, 0) #binary values 1/0 (selected/not selected)
lasso_matrix <- lasso_matrix[2:nrow(lasso_matrix),]
stable_variables <- as.data.frame(rowSums(lasso_matrix)) #data frame with count of how many times each variable is selected for a model
stable_variables$protein <- rownames(stable_variables) #column of variable names
colnames(stable_variables)[1] <- "times"
stable_variables <- stable_variables[!is.na(stable_variables$times),]  #remove NAs
stable_variables <- stable_variables[stable_variables$times != 0,] #remove all not selected variables
stable_variables <- stable_variables[order(-stable_variables$times),] #ordering by times 

###plotting 50 stable variables
(ggplot(stable_variables[1:50,], aes(x=reorder(as.factor(protein),-abs(times),mean), 
                                    y=times, col =reorder(as.factor(protein),-abs(times),mean), 
                                    fill =reorder(as.factor(protein),-abs(times),mean))) 
  + geom_col(show.legend = FALSE, alpha = 0.6) + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                                       axis.text=element_text(size=10)) 
  + xlab("Protein ID") + ylab("Times selected")
  + ggtitle("Stable variables obtained by Lasso regression"))
ggsave('3. Feature selection/BC_LassoProteins_lassoRegression.pdf')

stab_var <- stable_variables$protein[1:50] 
stab_var.ind <- which(colnames(data) %in% stab_var)

###boxplots for first 50 proteins obtained by lasso regression
plot_protein<- 1:length(stab_var[1:50])
plot_list <- list()
for (i in seq_along(plot_protein)){ 
  plot_list[[length(plot_list)+1]]  <- ggplot(data, aes_string("PAM50.mRNA", stab_var[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + 
    geom_boxplot(alpha=0.3) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
plot_grid(plotlist = plot_list,ncol = 5, nrow = 10)
ggsave("3. Feature selection/Lasso_proteins_boxplot.pdf", width = 70, height = 100, units = "cm")

###violin plots for first 50 proteins obtained by lasso regression
plot_protein<- 1:length(stab_var[1:50])
plot_list <- list()
for (i in seq_along(plot_protein)){ 
  plot_list[[length(plot_list)+1]]  <- ggplot(data, aes_string("PAM50.mRNA", stab_var[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + 
    geom_violin(alpha=0.3) + geom_boxplot(width=0.1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
plot_grid(plotlist = plot_list,ncol = 5, nrow = 10)
ggsave("3. Feature selection/Lasso_proteins_violinplot.pdf", width = 70, height = 100, units = "cm")

###create csv file with information about 50 proteins obtained by lasso regression 
library(biomaRt) #biomaRt database
listEnsembl()
biomart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
searchAttributes(mart = biomart) #search attributes
getBM(attributes=c('external_gene_name','refseq_peptide','entrezgene_description'),filters='refseq_peptide',
      values=stab_var, mart=biomart)->stab_var.biomart # biomaRt query
stable_variables.50 <- stable_variables[1:50,]
filter(stable_variables.50, stable_variables.50$protein %in% stab_var.biomart$refseq_peptide) -> time
colnames(time) <- c('times_lassoRegression', 'refseq_peptide')
stab_var.biomart <- merge(stab_var.biomart, time, by.x = 'refseq_peptide')
stab_var.biomart[order(stab_var.biomart$times_lassoRegression, decreasing = TRUE), ] -> stab_var.biomart
write.table(stab_var.biomart, "3. Feature selection/Lasso_proteins.csv", sep = ',', row.names = F)

#---------------------------------------------------------------------------------

## Creation the train and test data ###

###test and train split index
set.seed(1)
samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)
train_data <- data[samp, c(21, stab_var.ind)]
summary(train_data)
test_data <- data[-samp, c(21, stab_var.ind)]
summary(test_data)

###test and train data barplot
train_data.plot <- train_data
train_data.plot$data <- "train_data"
test_data.plot <- test_data
test_data.plot$data <- "test_data"
train_data.plot <- t(train_data.plot)
test_data.plot <- t(test_data.plot)
data.plot <- cbind.data.frame(train_data.plot,test_data.plot)
data.plot <- t(data.plot)
data.plot <- as.data.frame(data.plot)
ggplot(data.plot, aes(x=PAM50.mRNA, fill=data)) + geom_bar() +  ggtitle("train_data (0.7) and test_data (0.3) for Lasso proteins")
ggsave('4. Models buildings for Lasso proteins/Lasso_train&testData.pdf')

#--------------------------------------------------------------------------------

## Models Building ##

###Support-Vector Machine
set.seed(1)
train_control <- trainControl(method="repeatedcv", 
                              number=3, 
                              repeats=10, 
                              savePredictions = TRUE, 
                              summaryFunction = multiClassSummary) #settng up train control for cross validation and calibration of hyperparameter
grid <- expand.grid(C = seq(0.000001,0.15,0.002)) #Tunegrid for different values of C
set.seed(1)
svm.lin.mod <- train(PAM50.mRNA ~ ., 
                     data=train_data, 
                     trControl=train_control, 
                     method="svmLinear", 
                     preProcess = c("center","scale"), 
                     tuneGrid =grid, 
                     tuneLength = 10) #training
svm.predicts <- predict(svm.lin.mod, 
                        newdata = test_data) #Creating predictions on test set
sink("4. Models buildings for Lasso proteins/Lasso_SVM_results.txt", append =TRUE) #write output to file
cat("SVM Classifier for Lasso proteins")
cat("\n")
confusionMatrix(svm.predicts, as.factor(test_data$PAM50.mRNA)) #confusion matrix
sink() #stop recording output

###Random Forest
set.seed(1)
random.forest.mod <- randomForest(as.factor(PAM50.mRNA) ~ ., 
                                  data=train_data, 
                                  ntree = 500, 
                                  mtry = 8, 
                                  importance = TRUE) #training
plot(random.forest.mod)
random.forest.predicts <- predict(random.forest.mod, newdata=test_data, type = "class") #prediction
sink("4. Models buildings for Lasso proteins/Lasso_RandomForest_results.txt", append =TRUE) #write output to file
cat("RandomForest Classifier for Lasso proteins")
cat("\n")
confusionMatrix(random.forest.predicts, as.factor(test_data$PAM50.mRNA)) #confusion matrix
sink() #stop recording output

###Neural Network Classifier
neural.network.mod <-neuralnet(as.factor(PAM50.mRNA) ~ .,          
                               data=train_data, 
                               hidden = c(30, 15, 8), 
                               linear.output = FALSE,
                               act.fct = "logistic",
                               algorithm = 'sag')
plot(neural.network.mod,
     intercept = F,
     rep = 'best',
     col.hidden = 'darkblue', 
     col.hidden.synapse = 'black', 
     show.weights = F, 
     fill = 'lightblue',
     col.intercept	= 'red',
     fontsize = 8, 
     dimension = 15)
neural.network.mod$result.matrix #errors for model
neural.network.predicts <- compute(neural.network.mod, test_data) #prediction
neural.network.res = neural.network.predicts$net.result #results of prediction
print(neural.network.res)
neural.network.res=data.frame("neural.network.res"=ifelse(max.col(neural.network.res[ ,1:4])==1, "Basal-like",
                              ifelse(max.col(neural.network.res[ ,1:4])==2, "HER2-enriched",
                              ifelse(max.col(neural.network.res[ ,1:4])==3, "Luminal A", "Luminal B"))))
sink("4. Models buildings for Lasso proteins/Lasso_NeuralNetworkClassifier_results.txt", append =TRUE) #write output to file
cat("Neural Network Classifier for Lasso proteins")
cat("\n")
confusionMatrix(as.factor(test_data$PAM50.mRNA), as.factor(neural.network.res$neural.network.res))
sink() #stop recoridng output

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# PAM50 Models #

## Data modification, train and test data creation

###identify PAM50 proteins
pam50 <- genes$RefSeqProteinID
pam50.ind <- which(colnames(data) %in% pam50 )
length(pam50.ind)

###dividing the data set into a training data and a test data - 70% to 30%
samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)
train_data_pam50 <- data[samp, c(21, pam50.ind)]
summary(train_data_pam50)
test_data_pam50 <- data[-samp, c(21, pam50.ind)]
summary(test_data_pam50)

###test and train data barplot
train_data_pam50.plot <- train_data_pam50
train_data.plot_pam50$data <- "train_data_pam50"
test_data_pam50.plot <- test_data_pam50
test_data_pam50.plot$data <- "test_data_pam50"
train_data_pam50.plot <- t(train_data_pam50.plot)
test_data_pam50.plot <- t(test_data_pam50.plot)
data_pam50.plot <- cbind.data.frame(train_data_pam50.plot,test_data_pam50.plot)
data_pam50.plot <- t(data_pam50.plot)
data_pam50.plot <- as.data.frame(data_pam50.plot)
ggplot(data.plot, aes(x=PAM50.mRNA, fill=data)) + geom_bar() +  ggtitle("train_data (0.7) and test_data (0.3) for PAM50 proteins")
ggsave('5. Models buildings for PAM50 proteins/PAM50_train&testData.pdf')

###boxplots for PAM50 proteins
data_pam50 <- rbind(train_data_pam50,test_data_pam50)
data_protein <- data_pam50[,2:length(data_pam50)]
plot_protein <- colnames(data_protein)
plot_list <- list()
for (i in seq_along(plot_protein)){ 
  plot_list[[length(plot_list)+1]]  <- ggplot(data_pam50, aes_string("PAM50.mRNA", plot_protein[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + 
    geom_boxplot(alpha=0.3) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
plot_grid(plotlist = plot_list,ncol = 5, nrow = 9)
ggsave("5. Models buildings for PAM50 proteins/PAM50_proteins_boxplot.pdf", width = 70, height = 100, units = "cm")

###violin plots for PAM50 proteins
plot_protein <- colnames(data_protein)
plot_list <- list()
for (i in seq_along(plot_protein)){ 
  plot_list[[length(plot_list)+1]]  <- ggplot(data, aes_string("PAM50.mRNA", plot_protein[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + 
    geom_violin(alpha=0.3) + geom_boxplot(width=0.1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
plot_grid(plotlist = plot_list,ncol = 5, nrow = 9)
ggsave("5. Models buildings for PAM50 proteins/PAM50_proteins_violinplot.pdf", width = 70, height = 100, units = "cm")

#--------------------------------------------------------------------------------

## Model Building

###Support-Vector Machine
set.seed(1)
train_control_pam50 <- trainControl(method="repeatedcv", 
                              number=3, 
                              repeats=10, 
                              savePredictions = TRUE, 
                              summaryFunction = multiClassSummary) #settng up train control for cross validation and calibration of hyperparameter
grid <- expand.grid(C = seq(0.000001,0.15,0.002)) #Tunegrid for different values of C
set.seed(1)
svm.lin.mod_pam50 <- train(PAM50.mRNA ~ ., 
                           data=train_data_pam50, 
                           trControl=train_control_pam50, 
                           method="svmLinear", 
                           preProcess = c("center","scale"), 
                           tuneGrid =grid, 
                           tuneLength = 10) #training
svm.predicts_pam50 <- predict(svm.lin.mod_pam50, 
                              newdata = test_data_pam50) #predictions on test set
sink("5. Models buildings for PAM50 proteins/PAM50_SVM_results.txt", append =TRUE) #write output to file
cat("SVM Classifier for PAM50 proteins")
cat("\n")
confusionMatrix(svm.predicts_pam50, as.factor(test_data_pam50$PAM50.mRNA)) # confusion matrix
sink() #stop recoridng output

###Random Forest
set.seed(1)
random.forest.mod_pam50 <- randomForest(as.factor(PAM50.mRNA) ~ ., 
                                        data=train_data_pam50, 
                                        ntree = 500, 
                                        mtry = 8, 
                                        importance = TRUE) #training
plot(random.forest.mod_pam50)
random.forest.predicts_pam50 <- predict(random.forest.mod_pam50, 
                                        newdata=test_data_pam50, 
                                        type = "class") #prediction
confusionMatrix(random.forest.predicts_pam50, as.factor(test_data_pam50$PAM50.mRNA)) #confusion matrix
sink("5. Models buildings for PAM50 proteins/PAM50_RandomForest_results.txt", append =TRUE) #write output to file
cat("Random Forest Classifier for PAM50 proteins")
cat("\n")
confusionMatrix(random.forest.predicts_pam50, as.factor(test_data_pam50$PAM50.mRNA)) #confusion matrix
sink() #stop recoridng output

###Neural Network Classifier
neural.network.mod_pam50 <-neuralnet(as.factor(PAM50.mRNA) ~ .,          
                                    data=train_data_pam50, 
                                    hidden = c(30, 15, 8), 
                                    linear.output = FALSE,
                                    act.fct = "logistic",
                                    algorithm = 'sag')
plot(neural.network.mod_pam50,
     intercept = F,
     rep = 'best',
     col.hidden = 'darkblue', 
     col.hidden.synapse = 'black', 
     show.weights = F, 
     fill = 'lightblue',
     col.intercept	= 'red',
     fontsize = 8, 
     dimension = 15)
neural.network.mod_pam50$result.matrix #errors for model
neural.network.predicts_pam50 <- compute(neural.network.mod_pam50, test_data_pam50) #prediction
neural.network.res_pam50 = neural.network.predicts_pam50$net.result #results of prediction
print(neural.network.res_pam50)
neural.network.res_pam50=data.frame("neural.network.res_pam50"=ifelse(max.col(neural.network.res_pam50[ ,1:4])==1, "Basal-like",
                                                          ifelse(max.col(neural.network.res_pam50[ ,1:4])==2, "HER2-enriched",
                                                                 ifelse(max.col(neural.network.res_pam50[ ,1:4])==3, "Luminal A", "Luminal B"))))
sink("5. Models buildings for PAM50 proteins/PAM50_NeuralNetworkClassifier_results.txt", append =TRUE) #write output to file
cat("Neural Network Classifier for PAM50 proteins")
cat("\n")
confusionMatrix(as.factor(test_data_pam50$PAM50.mRNA), as.factor(neural.network.res_pam50$neural.network.res_pam50))
sink() #stop recoridng output

