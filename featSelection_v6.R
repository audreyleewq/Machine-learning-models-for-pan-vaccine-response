# HIPC2 Signatures Project
# Leave-one-vaccine-out CV antibody response prediction pipeline with feature selection
# Added logistic and lasso option for feature selection and model prediction functions
# 5/29/2020 - added glmnet function instead of just calling cv.glmnet. Added code to get only upGenes across all pt as features


library(RColorBrewer)
library(ggthemes)
library(stats)
library(tidyverse)
library(dplyr)
library(reshape2)
library(qvalue)
library(genefilter)
library(pheatmap)
library(Biobase)
library(DescTools)
library(matrixStats)
library(corrplot)
library(org.Hs.eg.db)
library(fgsea)
library(data.table)
library(glmnet)
library(pROC)
library(leaps)
library(caret)

rm(list=ls())
# #Audrey paths
setwd("/Users/audrey/Desktop/hipc/")
#Load master color set
load("adj_path_vt_colors.rda")
#Load GE
eset <- readRDS("2020_05_28_IS2_eset_withResponse_norm_young.rds")
#Load BTMs
load('BTM_for_GSEA_20131008_geneID.RData')
#Thomas paths
# setwd("/Users/tlhagan/Documents/Stanford/HIPC 2/Signatures Project/Analysis/Virtual Study/")
# #Load master color set
# load("adj_path_vt_colors.rda")
# #Load GE
# load("CURRENT 2020_02_18/2020_02_18_IS2_withResponse_normCorr.rda")
# #Load BTMs
# load('~/Documents/Emory/BTM_files/BTM_for_GSEA_20131008_geneID.RData')

#Add subject colnames
colnames(eset)=colnames(exprs(eset))
#Cutoff/input options
p_cutoff=0.05
fdr_cutoff=0.05
NES_cutoff=2
corrMethod='spearman'
BTM = TRUE
ab_response_col='MFC'
responder_col='MFC_p30'
feat_select_opt='binary' 
model_select='lasso' #logistic, lasso
n_genes = 150 #how many genes to start with for feature selection
n_feat_select = 25  #maximum number of features - must be smaller than n_genes
increaseThreshold=0
nRounds=5  #stop feature selection after nRounds of training auc not improving
add_age_gender=FALSE
combining=FALSE  #Should we combine to 'other' vaccines


########### FUNCTIONS ############

#Function to pick best kth model for each k predictor based on training auc. 
#Take as input trainData as eset, totalFeatures = n_feat_select (max k), modelSelect = 'logistic' or 'lasso'
#Returns a list with 2 sublists - bestFeatures[[1]] each member contains a vector of best feature indices
#bestFeatures[[2]] each member contains the training AUC for the selected features model
select_predictors <- function(trainData, totalFeatures, modelSelect, nRounds=5, increaseThreshold=0) {
  set.seed(123)
  #Initialize bestFeatures list 
  bestFeatures <- vector("list", 2) %>% set_names(c("feature_indices", "train_auc")) 
  for(index in 1:length(bestFeatures)) {bestFeatures[[index]] <- vector("list", totalFeatures)}
  bestFeatures[[1]][[1]] <- 1  #start with first predictor (no p=1 model because lasso takes min. p=2)
  #Find best features for each p=k model  
  modelTrainScore <- vector("list", totalFeatures)
  modelTrainScore[[1]] <- 0
  counter=0  #counts number of times training auc did not improve
  for(n in 2:totalFeatures) {  #start with 2 because we start with p=2 model
    print(paste0('k=',n))
    bestTrainScore <- 0
    #Try every predictor for k number of predictors
    for(i in 1:nrow(trainData)) {
      if(i %in% bestFeatures[[1]][[n-1]]) next  #skip if feature is already selected
      print(i)
      currentFeats <- c(bestFeatures[[1]][[n-1]], i)  #add 1 feature
      trainDataSubset <- trainData[currentFeats,]
      #Perform feature selection using chosen model method - lasso or logistic
      if(modelSelect == 'lasso') {
        grid =10^ seq(10,-2,length=100)  
        model = glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha=1, lambda=grid, weights=trainDataSubset$weight, family='binomial')
        cv_out = cv.glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha = 1) #find the best lambda
        best_lam = cv_out$lambda.min
        # Predict training auc. "response" gives the probability 
        predProbs <- predict(model, s=best_lam, newx=as(t(exprs(trainDataSubset)),"dgCMatrix"), type="response")
      } 
      if(modelSelect == 'logistic') {
        model <- glmnet(t(exprs(trainDataSubset)),trainDataSubset$response_cat, weights=trainDataSubset$weight, 
                        family='binomial', lambda = 0)
        predProbs <- predict(model, s=0, newx=as(t(exprs(trainDataSubset)),"dgCMatrix"), type="response")
      }
      predClass = ifelse(predProbs>0.5, 1, 0)  
      currTrainScore <- auc(trainDataSubset$response_cat, as.vector(predClass))
      if(currTrainScore > bestTrainScore) {
        bestTrainScore <- currTrainScore
        bestFeatures[[1]][[n]] <- currentFeats
        bestFeatures[[2]][[n]] <- currTrainScore
      }
    }
    #Stop feature selection when training auc has not improved for nRounds
    modelTrainScore[[n]] <- bestTrainScore
    if(bestTrainScore > (modelTrainScore[[n-1]]+increaseThreshold)) {counter=0} else {counter=counter+1}
    if(counter==nRounds) {
      bestFeatures <- lapply(bestFeatures, function(x) x[1:(n-counter)])
      break
    }
  }
  return(bestFeatures)
}

# Function to run test over different training CV sets. 
# Take as input training data eset, number of fold CV, bestFeatures vector, modelSelect='logistic' or 'lasso'
# Output: a list with nFolds validation scores for each feature subset
get_test_varability <- function(train_input, nFolds, bestFeatures, test_input, modelSelect) {
  set.seed(123)
  folds = createFolds(train_input[[responder_col]], k=nFolds, list=TRUE, returnTrain=FALSE)
  test_scores <- c()
  for (k in 1:nFolds) {  #For each CV set
    #Split into training/validation
    train_cv_input=train_input[,-folds[[k]]]
    val_cv_input=train_input[,folds[[k]]]
    #Get a list of test scores each calculated from a subset of training data model 
    test_result <- model_test_result(train_cv_input, test_input, bestFeatures, modelSelect)
    test_scores <- c(test_scores, test_result)
  }
  return(test_scores)  
}

#Function to train and test model based on selected features
#Take as input train and test data eset, selected feature indices, modelSelect='logistic' or 'lasso'
#Output: An AUC score for the test data
model_test_result <- function(trainData, testData, featureInd, modelSelect) {
  set.seed(123)
  trainDataSubset <- trainData[featureInd,]
  testDataSubset <- testData[featureInd,]
  # Train data - should we use type.measure="auc"?
  if(modelSelect == 'lasso') {
    grid =10^ seq(10,-2,length=100)  
    model = glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha=1, lambda=grid, weights=trainDataSubset$weight, family='binomial')
    cv_out = cv.glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha = 1) #find the best lambda
    best_lam = cv_out$lambda.min
    # Predict test data
    predProbs <- predict(model, s=best_lam, newx=as(t(exprs(testDataSubset)),"dgCMatrix"), type="response")
  }
  if(modelSelect == 'logistic') {
    model <- glmnet(t(exprs(trainDataSubset)),trainDataSubset$response_cat, weights=trainDataSubset$weight, 
                    family='binomial', lambda = 0)
    predProbs <- predict(model, s=0, newx=as(t(exprs(testDataSubset)),"dgCMatrix"), type="response")
  }
  predClass = ifelse(predProbs>0.5, 1, 0)
  testScore <- auc(testDataSubset$response_cat, as.vector(predClass))
  return(testScore)
}

##For getting test prediction accuracy across CVs
get_test_accuracy <- function(train_input, nFolds, bestFeatures, test_input, modelSelect) {
  set.seed(123)
  folds = createFolds(train_input[[responder_col]], k=nFolds, list=TRUE, returnTrain=FALSE)
  test_scores <- c()
  for (k in 1:nFolds) {  #For each CV set
    #Split into training/validation
    train_cv_input=train_input[,-folds[[k]]]
    val_cv_input=train_input[,folds[[k]]]
    #Get a list of test scores each calculated from a subset of training data model 
    test_result <- model_prediction_accuracy(train_cv_input, test_input, bestFeatures, modelSelect)
    test_scores <- c(test_scores, test_result)
  }
  return(test_scores)  
}

#Get prediction accuracy 
model_prediction_accuracy <- function(trainData, testData, featureInd, modelSelect) {
  set.seed(123)
  trainDataSubset <- trainData[featureInd,]
  testDataSubset <- testData[featureInd,]
  # Train data - should we use type.measure="auc"?
  if(modelSelect == 'lasso') {
    grid =10^ seq(10,-2,length=100)  
    model = glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha=1, lambda=grid, weights=trainDataSubset$weight, family='binomial')
    cv_out = cv.glmnet(t(exprs(trainDataSubset)), trainDataSubset$response_cat, alpha = 1) #find the best lambda
    best_lam = cv_out$lambda.min
    # Predict test data
    predProbs <- predict(model, s=best_lam, newx=as(t(exprs(testDataSubset)),"dgCMatrix"), type="response")
  }
  if(modelSelect == 'logistic') {
    model <- glmnet(t(exprs(trainDataSubset)),trainDataSubset$response_cat, weights=trainDataSubset$weight, 
                    family='binomial', lambda = 0)
    predProbs <- predict(model, s=0, newx=as(t(exprs(testDataSubset)),"dgCMatrix"), type="response")
  }
  predClass = ifelse(predProbs>0.5, 1, 0)
  testScore <- length(which(testDataSubset$response_cat == as.vector(predClass)))/length(testDataSubset$response_cat)
  return(testScore)
}


############### TIDY DATASET ###################

#Remove 'viral vector' from recombinant protein labels
eset$vaccine_type=recode(eset$vaccine_type,"Recombinant protein+viral vector"="Recombinant protein")
#Reclassify live-attenuated vaccines to be adjuvanted
levels(eset$adjuvant)=c(levels(eset$adjuvant), 'Live attenuated')
ind=grepl('Live attenuated', eset$vaccine_type)
eset$adjuvant[ind]='Live attenuated'
#Create combined vaccine type/pathogen column
eset$pathogen_type=paste(eset$pathogen," (",eset$vaccine_type,")")

#Should we combine vaccines to 'others'?
if(combining) {
  othersInd <- which(eset$pathogen != "Influenza")
  eset$pathogen_type[othersInd] <- "Others"
}

#Prune list
#Remove subjects with Inf/NA MFC measurements
eset=eset[,!is.na(eset$MFC)]  #study1260 meningococcus got no MFC so remove NA
# eset=eset[,-which(eset$MFC == Inf)]  #exclude study1289 which has MFC = Inf
#Right now SDY1260 has MFC values but they are incorrect-need to temporarily remove
# eset=eset[,-which(eset$study2=='SDY1260')]

#Find samples from timepoints of interest
tp_int=c(0,3,7) #Only using Days 3 and 7 for now
#tp_int=unique(eset$study_time_collected[which(eset$study_time_collected>=0)]) #Alternate: use all timepoints (>=0)
ind=lapply(tp_int, function(x) which(eset$study_time_collected==x))
#Retain only samples from timepoints of interest
eset=eset[,Reduce(union,ind)]
#Recompute timepoint indices after removing extraneous timepoints
ind=lapply(tp_int, function(x) which(eset$study_time_collected==x))
#Remove samples from a single study with fewer than sample_cutoff samples at any timepoint
sample_cutoff = 3
matrix_uni_tp=lapply(ind,function(x) unique(eset[,x]$matrix))
matrix_ind=lapply(1:length(ind),function(x)
  lapply(1:length(matrix_uni_tp[[x]]), function(y) which(matrix_uni_tp[[x]][[y]]==eset[,ind[[x]]]$matrix)))
ind_cut_all=vector()
for (i in 1:length(matrix_ind)) {
  ind_cut=which(sapply(matrix_ind[[i]],length)<sample_cutoff)
  if (is_empty(ind_cut)==FALSE) {
    #edited this line because new VS has >1 matrix in ind cut - SDY1119 has only 1 sample in the new VS
    # ind_total = sapply(ind_cut, function(x) ind[[i]][matrix_ind[[i]][[x]]])
    # ind_cut_all=c(ind_cut_all, ind_total)  
    ind_cut_all=c(ind_cut_all,ind[[i]][matrix_ind[[i]][[ind_cut]]])
  }
}
if (is_empty(ind_cut_all)==FALSE) {
  eset=eset[,-ind_cut_all]
}

#Recompute timepoint indices after removing samples
tp_int=unique(eset$study_time_collected[which(eset$study_time_collected>=0)])
ind=lapply(tp_int, function(x) which(eset$study_time_collected==x))

#Remove genes with NA
eset=eset[complete.cases(exprs(eset)),]

#Create combined SDY/pathogen/vaccine type column
eset$SDY_pathogen_type=paste(eset$study,eset$pathogen_type)
#Create unique list of studies
matrix_uni=unique(eset$matrix)

#Create master color list for pathogen/vaccine type combination
cols=colorRampPalette(brewer.pal(length(unique(eset$pathogen_type)), "Set3"))
mycolors=cols(length(unique(eset$pathogen_type)))
names(mycolors)=unique(eset$pathogen_type)

#If BTM inputs are selected, collapse to BTM level
if (BTM) {
  #Convert genes to entrez ID
  entrez_gene_lookup=select(org.Hs.eg.db, rownames(eset), c("ENTREZID"), "ALIAS")
  #Remove genes with multiple entrez IDs
  entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,2]) | duplicated(entrez_gene_lookup[,2], fromLast=TRUE)),]
  #Remove entrez IDs with multiple genes
  entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,1]) | duplicated(entrez_gene_lookup[,1], fromLast=TRUE)),]
  entrez_list=as.character(entrez_gene_lookup[match(rownames(eset),entrez_gene_lookup[,1]),2])  #get entrezID of genes that are in eset
  #Remove BTMs which have no matching genes in dataset
  BTM_list=BTM_list[!sapply(BTM_list, function(x) is_empty(intersect(x,entrez_list)))]
  #Collapse - match entrezID in BTM_list to entrezID of genes (arrange according to rownames of eset) - find arithmetic mean for each sample
  exp_BTM=do.call(rbind, lapply(BTM_list, function(x) colMeans(exprs(eset[na.omit(match(x,entrez_list)),]),na.rm=TRUE)))
  #Create BTM eset (maybe not proper approach but it works)
  eset_BTM=eset[1:nrow(exp_BTM),]
  rownames(eset_BTM)=rownames(exp_BTM)
  exprs(eset_BTM)=exp_BTM
  #Replace eset with BTM eset
  eset=eset_BTM
  rm(eset_BTM)
}

#Compute D0 normalized FC - for individual genes/BTMs - exp_FC has 3 large matrix with D1, D3, D7 timepoints
ind_D0=which(0==eset$study_time_collected)
common=lapply(2:length(ind),function(x) intersect(eset$participant_id[ind[[x]]],eset$participant_id[ind_D0]))  #get participant id with timepoint and D0 data
ia=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind[[x]]]))) #included participant IDs match with participant IDs at diff timepoints in eset
ib=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind_D0])))  #participant IDs match with D0 to get indices in ind_D0
exp_FC=lapply(2:length(ind),function(x) eset[,ind[[x]][ia[[x-1]]]])  #get expression data of the 3 timepoints
exp_FC=lapply(2:length(ind),function(x) {exprs(exp_FC[[x-1]])=exprs(exp_FC[[x-1]])-exprs(eset[,ind_D0[ib[[x-1]]]]); exp_FC[[x-1]]})

#For each study, average expression across all subjects
#Study matrices for each timepoint
matrix_uni_tp=lapply(exp_FC,function(x) x$matrix[!duplicated(x$matrix)])  #study matrix for each timepoint
#Store study metadata
matrix_uni_tp_metaData=lapply(exp_FC,function(x) pData(x)[!duplicated(x$matrix),] %>%
                                rownames_to_column() %>%
                                dplyr::select(matrix, vaccine, vaccine_type, pathogen, adjuvant, pathogen_type) %>%
                                column_to_rownames(var = "matrix"))
#Indices for each study
matrix_ind=lapply(1:length(exp_FC),
                  function(x) lapply(1:length(matrix_uni_tp[[x]]), function(y) which(matrix_uni_tp[[x]][[y]]==exp_FC[[x]]$matrix)))
#High responder indices within each study
matrix_HR_ind=lapply(1:length(exp_FC),
                     function(x) lapply(1:length(matrix_uni_tp[[x]]),
                                        function(y) which(pData(exp_FC[[x]][,matrix_ind[[x]][[y]]])[responder_col]=='highResponder')))
matrix_LR_ind=lapply(1:length(exp_FC),
                     function(x) lapply(1:length(matrix_uni_tp[[x]]),
                                        function(y) which(pData(exp_FC[[x]][,matrix_ind[[x]][[y]]])[responder_col]=='lowResponder')))
#Study size
exp_FC_n=lapply(1:length(exp_FC),
                function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) length(matrix_ind[[x]][[y]])))

#Average FC
exp_FC_mean=lapply(1:length(exp_FC),
                   function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) rowMeans(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]))))

#Correlation with Ab response
# exp_FC_corr_stat=
#   lapply(1:length(exp_FC),  #for each timepoint
#          function(x) sapply(1:length(matrix_uni_tp[[x]]),  #for each study
#                             function(y) apply(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  #for each gene
#                                               function(z) cor.test(z, pData(exp_FC[[x]][,matrix_ind[[x]][[y]]])[[ab_response_col]], method = corrMethod)$statistic)))

#T-test between high/low responder
#If there are less than 2 HR/LR in a given study, return NA for all t-tests
exp_FC_t_stat=
  lapply(1:length(exp_FC),
         function(x) sapply(1:length(matrix_uni_tp[[x]]),
                            function(y) if(length(matrix_HR_ind[[x]][[y]])>=2 & length(matrix_LR_ind[[x]][[y]])>=2) apply(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]), 1, function(z) t.test(z[matrix_HR_ind[[x]][[y]]], z[matrix_LR_ind[[x]][[y]]])$statistic) else rep(NA, nrow(exp_FC[[x]]))))

#Find unique pathogen/vaccine types by timepoint
pt_uni_tp=lapply(matrix_uni_tp_metaData,function(x) x$pathogen_type[!duplicated(x$pathogen_type)])
pt_metaData=lapply(1:length(matrix_uni_tp_metaData),
                   function(x) dplyr::select(
                     matrix_uni_tp_metaData[[x]][!duplicated(matrix_uni_tp_metaData[[x]]$pathogen_type),],
                     vaccine_type, pathogen, adjuvant))

#Sample indices per vaccine
pt_ind=lapply(1:length(exp_FC),
              function(x) sapply(1:length(pt_uni_tp[[x]]),
                                 function(y) which(pt_uni_tp[[x]][[y]]==matrix_uni_tp_metaData[[x]]$pathogen_type)))

#Merge mean FC, correlation/t-test test statistics by vaccine (weighted average by study size)
pt_meanFC=
  lapply(1:length(exp_FC),
         function(x) sapply(1:length(pt_uni_tp[[x]]),
                            function(y) if(length(pt_ind[[x]][[y]])>1) rowWeightedMeans(exp_FC_mean[[x]][,pt_ind[[x]][[y]]], w=exp_FC_n[[x]][pt_ind[[x]][[y]]], na.rm=TRUE) else exp_FC_mean[[x]][,pt_ind[[x]][[y]]]))

# pt_corr_stat=
#   lapply(1:length(exp_FC),
#          function(x) sapply(1:length(pt_uni_tp[[x]]),
#                             function(y) if(length(pt_ind[[x]][[y]])>1) rowWeightedMeans(exp_FC_corr_stat[[x]][,pt_ind[[x]][[y]]], w=exp_FC_n[[x]][pt_ind[[x]][[y]]], na.rm=TRUE) else exp_FC_corr_stat[[x]][,pt_ind[[x]][[y]]]))

pt_t_stat=
  lapply(1:length(exp_FC),
         function(x) sapply(1:length(pt_uni_tp[[x]]),
                            function(y) if(length(pt_ind[[x]][[y]])>1) rowWeightedMeans(exp_FC_t_stat[[x]][,pt_ind[[x]][[y]]], w=exp_FC_n[[x]][pt_ind[[x]][[y]]], na.rm=TRUE) else exp_FC_t_stat[[x]][,pt_ind[[x]][[y]]]))


################ PREDICTION ######################
set.seed(123)
features_selected <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_auc_all <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pt_top_val_scores <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_auc_average <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_auc_sd <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_accuracy_all <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_accuracy_average <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
test_accuracy_sd <- vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
bestFeatures <- list()

# For each timepoint, take training data
for (i in 1:length(exp_FC)) {
  #Set names for lists to keep
  features_selected[[i]] <- vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pt_top_val_scores[[i]] <- vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  
  #For each pathogen/type, separate that vaccine as test set and use all others as training
  for (j in 1:length(pt_uni_tp[[i]])) {    ##change pt_uni_tp to combine vaccines to 'others'
    #Indices of samples from test vaccine
    pt_test_ind=which(pt_uni_tp[[i]][[j]]==exp_FC[[i]]$pathogen_type)
    #Indices of samples from training vaccines
    pt_train_ind=setdiff(1:ncol(exp_FC[[i]]),pt_test_ind)
    # Subset eset for training and testing data
    train_input=exp_FC[[i]][,pt_train_ind]
    test_input=exp_FC[[i]][,pt_test_ind]
    
    # Clean up training data
    if (add_age_gender==TRUE) {
      train_input=ExpressionSet(rbind(exprs(train_input),Age=as.numeric(train_input$age_imputed),Gender=as.numeric(train_input$gender_imputed=='Male')),
                                phenoData=train_input@phenoData)
      test_input=ExpressionSet(rbind(exprs(test_input),Age=as.numeric(test_input$age_imputed),Gender=as.numeric(test_input$gender_imputed=='Male')),
                               phenoData=test_input@phenoData)
    }
    #Remove moderate responders
    train_input=train_input[,!train_input[[responder_col]]=='moderateResponder']
    test_input=test_input[,!test_input[[responder_col]]=='moderateResponder']
    #Add numeric output
    train_input$response_cat=as.numeric(train_input[[responder_col]]=='highResponder')
    test_input$response_cat=as.numeric(test_input[[responder_col]]=='highResponder')
    
    #Order features based on train input - code is made to only binary output for now
    if (feat_select_opt=='binary') { #If using binary HR/LR comparison
      #Perform unweighted average of HR/LR t-test statistic within training vaccines
      pt_train_mean_stat=rowMeans(pt_t_stat[[i]][,-j], na.rm=TRUE)   ##change pt_t_stat for combining vaccines
      #Order train input based on top features (mean t-statistic for the pathogen type in training set)
      orderedFeatInd <- order(-abs(pt_train_mean_stat))
      train_input <- train_input[orderedFeatInd,]
      # Subset training set to include number of top genes to start with
      if(!BTM) {train_input <- train_input[1:n_genes,]}
      
      #get genes that are upregulated in all the studies
      # fcSign <- apply(sign(pt_meanFC[[i]][,-j]), 1, function(x) all(x == 1))
      # orderedFeatInd <- which(fcSign == TRUE)
      # train_input <- train_input[orderedFeatInd,]
      
      # #Try using TIV genes for prediction
      # tiv_genes_list <-readRDS("tiv_genes_list_lasso.rds")
      # orderedFeatInd <- which(rownames(train_input) %in% tiv_genes_list[[i]])
      # train_input <- train_input[orderedFeatInd,]
    }

    #To remove vaccine bias, weight each training sample by inverse of vaccine sample size
    train_pt_weight <- 1/sapply(pt_uni_tp[[i]][-j], function(x) length(which(train_input$pathogen_type==x)))
    train_input$weight <- train_pt_weight[match(train_input$pathogen_type,names(train_pt_weight))]
    
    #Perform feature selection and get best p=k model
    # bestFeatures <- select_predictors(train_input, n_feat_select, modelSelect=model_select, nRounds=nRounds, increaseThreshold=increaseThreshold)  #a list of feature indices for each p=k model
    # topFeatures <- unlist(bestFeatures[[1]][length(bestFeatures[[1]])])
    # features_selected[[i]][[j]] <- topFeatures  #store top features as reference
    
    # If dont want to perform feature selection
    topFeatures <- rownames(train_input)
    features_selected[[i]][[j]] <- topFeatures  
    
    #Train and test using top selected features
    ##If there is only 1 response level, can't calculate auc, calculate test prediction accuracy
    if(length(unique(test_input$response_cat)) < 2) {
      test_auc_all[[i]][[j]] <- NA
      next
    }
    ##Calculate both auc and prediction accuracy for all 
    testAUC <- get_test_varability(train_input, nFolds=10, topFeatures, test_input, modelSelect=model_select)
    test_auc_all[[i]][[j]] <- testAUC  #store test performance (over nFold CV)
    test_auc_average[[i]][j]=mean(testAUC)
    test_auc_sd[[i]][j]=sd(testAUC)
    testAccuracy <- get_test_accuracy(train_input, nFolds=10, topFeatures, test_input, modelSelect=model_select)
    test_accuracy_all[[i]][[j]] <- testAccuracy
    test_accuracy_average[[i]][[j]]=mean(testAccuracy)
    test_accuracy_sd[[i]][j]=sd(testAccuracy)
    
    # #try not performing CV
    # test_auc_average[[i]][[j]] <- model_test_result(train_input, test_input, topFeatures, model_select)
    # test_accuracy_average[[i]][[j]] <- model_prediction_accuracy(train_input, test_input, topFeatures, model_select)
  }  
}

test_auc_average <- lapply(1:length(test_auc_average), function(x) test_auc_average[[x]] %>% set_names(pt_uni_tp[[x]]))
test_auc_sd <- lapply(1:length(test_auc_sd), function(x) test_auc_sd[[x]] %>% set_names(pt_uni_tp[[x]]))

# #Store the features selected for each vaccine that were used to run the model
sink("BTMs_pt_features.csv")
print(features_selected)
sink()

# #Plot training auc over each p=k model
# unlist(bestFeatures[[2]])
# bestFeatures[[2]][[1]]=NA
# dfplot <- data.frame("nPredictors"=c(1:length(bestFeatures[[2]])), "trainingAUC"=unlist(bestFeatures[[2]]))
# ggplot(aes(x=nPredictors, y=trainingAUC), data=dfplot) + geom_point() + geom_smooth()


################# PLOTTING #####################
for (i in 1:length(exp_FC)) {
  df=rbind('AUC.GE'=test_auc_average[[i]],'Std.GE'=test_auc_sd[[i]])
  df=data.frame(t(df))
  df$Vaccine=factor(rownames(df), levels=rownames(df))
  #Tidy table (possibly messy, could be cleaned)
  df=df %>% gather(key='AUC.Type',value='AUC',AUC.GE) %>%
    gather(key='Std.Type',value='Std',Std.GE) %>% filter(substring(AUC.Type, 4) == substring(Std.Type, 4))
  ggplot(df,aes(x=Vaccine,y=AUC,ymin=AUC-Std,ymax=AUC+Std,fill=AUC.Type))+
    geom_bar(stat="identity",width=0.7,position=position_dodge(0.7))+
    geom_errorbar(width=0.4,position=position_dodge(0.7), color="black")+
    scale_fill_manual(values=c('darkgrey','lightgrey'))+
    coord_cartesian(ylim=c(0,1))+
    #coord_flip()+
    xlab("")+
    ylab("AUC")+
    theme_bw()+
    theme_set(theme_bw() + theme(legend.key=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)))
  ggsave(paste('Test_AUC_Day',tp_int[i+1], '_day_BTM_Pt_lasso_wCV', '.pdf', sep=''), width = 7, height = 8)
}

#Get list of input features used for the model
# df_d3 <- data.frame(lapply(features_selected[[1]], "length<-", max(lengths(features_selected[[1]]))))
# df_d7 <- data.frame(lapply(features_selected[[2]], "length<-", max(lengths(features_selected[[2]]))))
# fwrite(x = df_d3, file = "pt_d3_feats.csv")
# fwrite(x = df_d7, file = "pt_d7_feats.csv")
# n_selected_feats <- lapply(1:length(features_selected), function(x) sapply(features_selected[[x]], length)) %>% set_names(tp_int[-1])
# nfeats_pt <- lapply(n_selected_feats, function(x) data.frame("nfeats" = x))
# 
# #Get sample size for each pathogen type
# nsamples_pt <- lapply(exp_FC, function(y) data.frame("n" = sapply(unique(eset$pathogen_type), 
#                                                    function(x) length(which(y$pathogen_type == x))))) %>% set_names(tp_int[-1])

