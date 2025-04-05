# HIPC2 Signatures Project
# Leave-one-vaccine-out CV antibody response prediction pipeline
# 2/12/20

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
eset <- readRDS("2020_04_24_IS2_eset_withResponse_norm_young.rds")
# eset <- readRDS("2020_04_13_IS2_eset_withResponse_norm_young.rds")
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
inputFeatures='genes'  #genes or BTM
ab_response_col='MFC'
responder_col='MFC_p30'
feat_select_opt='binary' 
model_select='logistic' #logistic, lasso
n_feat_select=100
add_age_gender=FALSE

#Remove 'viral vector' from recombinant protein labels
eset$vaccine_type=recode(eset$vaccine_type,"Recombinant protein+viral vector"="Recombinant protein")
#Reclassify live-attenuated vaccines to be adjuvanted
levels(eset$adjuvant)=c(levels(eset$adjuvant), 'Live attenuated')
ind=grepl('Live attenuated', eset$vaccine_type)
eset$adjuvant[ind]='Live attenuated'
#Create combined vaccine type/pathogen column
eset$pathogen_type=paste(eset$pathogen," (",eset$vaccine_type,")")
#Prune list
#Remove subjects with Inf/NA MFC measurements
eset=eset[,!is.na(eset$MFC)]  #study1260 meningococcus got no MFC so remove NA
# eset=eset[,-(which(eset$MFC == Inf))]  #exclude study1289 which has MFC = Inf - study already excluded in new VS
#Right now SDY1260 has MFC values but they are incorrect-need to temporarily remove - SDY1260 is fixed in new VS
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
sample_cutoff=3
matrix_uni_tp = lapply(ind, function(x) unique(eset[,x]$matrix))
matrix_ind = lapply(1:length(ind), function(x) lapply(1:length(matrix_uni_tp[[x]]), 
                                                      function(y) which(matrix_uni_tp[[x]][[y]]==eset[,ind[[x]]]$matrix)))
ind_cut_all=vector()
for (i in 1:length(matrix_ind)) {
  ind_cut=which(sapply(matrix_ind[[i]],length)<sample_cutoff)
  if (is_empty(ind_cut)==FALSE) {
    #edited this line because new VS has >1 matrix in ind cut - SDY1119 has only 1 sample in the new VS
    ind_total = sapply(ind_cut, function(x) ind[[i]][matrix_ind[[i]][[x]]])
    ind_cut_all=c(ind_cut_all, ind_total)  
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
if (inputFeatures=="BTM") {
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


#Perform leave-one-vaccine-out CV of Ab response prediction
selected_feat=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
selected_feat_ind=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
model=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_train=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_val=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_test=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_train=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_val=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_test=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_train_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_val_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_test_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_train=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_val=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_test=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_train_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_val_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_test_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
model_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_train_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_val_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
pred_test_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_train_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_val_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_test_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_train_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_val_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
accuracy_test_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_train_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_val_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_test_ag=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_train_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_val_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
auc_test_ag_std=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
best_lam_list=vector("list",length(exp_FC)) %>% set_names(tp_int[-1])
#For each timepoint
for (i in 1:length(exp_FC)) {
  selected_feat[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  selected_feat_ind[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  model[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_train[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_val[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_test[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_train[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_val[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_test[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_train_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_val_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_test_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_train[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_val[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_test[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_train_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_val_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_test_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  model_ag[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_train_ag[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_val_ag[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  pred_test_ag[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_train_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_val_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_test_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_train_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_val_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  accuracy_test_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_train_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_val_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_test_ag[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_train_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_val_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  auc_test_ag_std[[i]]=numeric(length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  best_lam_list[[i]]=vector("list",length(pt_uni_tp[[i]])) %>% set_names(pt_uni_tp[[i]])
  #For each pathogen/type, separate that vaccine as test set and use all others as training
  for (j in 1:length(pt_uni_tp[[i]])) {
    #Indices of samples from test vaccine
    pt_test_ind=which(pt_uni_tp[[i]][[j]]==exp_FC[[i]]$pathogen_type)
    #Indices of samples from training vaccines
    pt_train_ind=setdiff(1:ncol(exp_FC[[i]]),pt_test_ind)
    #Use training data to select features
    #Sort by abs value of test statistic of choice, select top features as inputs
    if (feat_select_opt=='binary') { #If using binary HR/LR comparison
      #Perform unweighted average of HR/LR t-test statistic within training vaccines
      pt_train_mean_stat=rowMeans(pt_t_stat[[i]][,-j], na.rm=TRUE)
      #Select top features
      selected_feat_ind[[i]][[j]]=order(-abs(pt_train_mean_stat))[1:n_feat_select]
      selected_feat[[i]][[j]]=names(pt_train_mean_stat)[order(-abs(pt_train_mean_stat))[1:n_feat_select]]
    } else { #If using correlation with antibody response
      #Perform unweighted average of correlation test statistic within training vaccines
      pt_train_mean_stat=rowMeans(pt_corr_stat[[i]][,-j], na.rm=TRUE)
      #Select top features
      selected_feat_ind[[i]][[j]]=order(-abs(pt_train_mean_stat))[1:n_feat_select]
      selected_feat[[i]][[j]]=names(pt_train_mean_stat)[selected_feat_ind[[i]][[j]]]
    }
    
    #Set training/testing inputs/outputs
    train_input=exp_FC[[i]][selected_feat_ind[[i]][[j]],pt_train_ind]
    test_input=exp_FC[[i]][selected_feat_ind[[i]][[j]],pt_test_ind]
    #Add age/gender as inputs if option is selected
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
    #Split into 10-fold cross-validation sets
    folds=createFolds(train_input[[responder_col]], k=10, list=TRUE, returnTrain=FALSE)
    model[[i]][[j]]=vector("list",length(folds))
    pred_train[[i]][[j]]=vector("list",length(folds))
    pred_test[[i]][[j]]=vector("list",length(folds))
    model_ag[[i]][[j]]=vector("list",length(folds))
    pred_train_ag[[i]][[j]]=vector("list",length(folds))
    pred_test_ag[[i]][[j]]=vector("list",length(folds))
    accuracy_train_cv=numeric(length(folds))
    accuracy_val_cv=numeric(length(folds))
    accuracy_test_cv=numeric(length(folds))
    accuracy_train_ag_cv=numeric(length(folds))
    accuracy_val_ag_cv=numeric(length(folds))
    accuracy_test_ag_cv=numeric(length(folds))
    auc_train_cv=numeric(length(folds))
    auc_val_cv=numeric(length(folds))
    auc_test_cv=numeric(length(folds))
    auc_train_ag_cv=numeric(length(folds))
    auc_val_ag_cv=numeric(length(folds))
    auc_test_ag_cv=numeric(length(folds))
    for (k in 1:length(folds)) {
      #Split into training/validation
      train_cv_input=train_input[,-folds[[k]]]
      val_cv_input=train_input[,folds[[k]]]
      
      #To remove vaccine bias, weight each training sample by inverse of vaccine sample size
      train_cv_pt_weight=1/sapply(pt_uni_tp[[i]][-j], function(x) length(which(train_cv_input$pathogen_type==x)))
      train_cv_input$weight=train_cv_pt_weight[match(train_cv_input$pathogen_type,names(train_cv_pt_weight))]
      
      ##Train classifier
      
      ##LASSO/Elastic Net 
      if(model_select == "lasso") {
        # grid =10^ seq(10,-2,length=100)  
        # model[[i]][[j]][[k]] = glmnet(t(exprs(train_cv_input)), train_cv_input$response_cat, alpha=1, lambda=grid, weights=train_cv_input$weight, family='binomial')
        # cv_out = cv.glmnet(t(exprs(train_cv_input)), train_cv_input$response_cat, alpha = 1) #find the best lambda
        # best_lam = cv_out$lambda.min
        model[[i]][[j]][[k]] = cv.glmnet(t(exprs(train_cv_input)), train_cv_input$response_cat, alpha = 1,
                                         weights=train_cv_input$weight, family='binomial') 
        best_lam = model[[i]][[j]][[k]]$lambda.min
        best_lam_list[[i]][[j]][[k]] = best_lam
        pred_train[[i]][[j]][[k]] = predict(model[[i]][[j]][[k]], s=best_lam, newx=as(t(exprs(train_cv_input)),"dgCMatrix"), type="response")
        pred_val[[i]][[j]][[k]] = predict(model[[i]][[j]][[k]], s=best_lam, newx=as(t(exprs(val_cv_input)),"dgCMatrix"), type="response")
        pred_test[[i]][[j]][[k]] = predict(model[[i]][[j]][[k]], s=best_lam, newx=as(t(exprs(test_input)),"dgCMatrix"), type="response")
        #Train using only age and gender
        # model_ag[[i]][[j]][[k]] = glmnet(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))), train_cv_input$response_cat, alpha=1, lambda=grid, weights=train_cv_input$weight, family='binomial')
        # cv_out = cv.glmnet(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))), train_cv_input$response_cat, alpha = 1) #find the best lambda
        # best_lam = cv_out$lambda.min
        model_ag[[i]][[j]][[k]] = cv.glmnet(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))), train_cv_input$response_cat, alpha = 1,
                                            weights=train_cv_input$weight, family='binomial') #find the best lambda
        best_lam = model_ag[[i]][[j]][[k]]$lambda.min
        pred_train_ag[[i]][[j]][[k]] = predict(model_ag[[i]][[j]][[k]], s=best_lam, newx=as(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))),"dgCMatrix"), type="response")
        pred_val_ag[[i]][[j]][[k]] = predict(model_ag[[i]][[j]][[k]], s=best_lam, newx=as(t(rbind(as.numeric(val_cv_input$age_imputed), as.numeric(val_cv_input$gender_imputed=='Male'))),"dgCMatrix"), type="response")
        pred_test_ag[[i]][[j]][[k]] = predict(model_ag[[i]][[j]][[k]], s=best_lam, newx=as(t(rbind(as.numeric(test_input$age_imputed), as.numeric(test_input$gender_imputed=='Male'))),"dgCMatrix"),type="response")
      }
      
      #Logistic Regression
      if(model_select == "logistic") {
        model[[i]][[j]][[k]]=glmnet(t(exprs(train_cv_input)),train_cv_input$response_cat, weights=train_cv_input$weight, family='binomial', lambda = 0)
        #Get model outputs for train and test data
        pred_train[[i]][[j]][[k]]=predict(model[[i]][[j]][[k]], newx=as(t(exprs(train_cv_input)),"dgCMatrix"), s=0)
        pred_val[[i]][[j]][[k]]=predict(model[[i]][[j]][[k]], newx=as(t(exprs(val_cv_input)),"dgCMatrix"), s=0)
        pred_test[[i]][[j]][[k]]=predict(model[[i]][[j]][[k]], newx=as(t(exprs(test_input)),"dgCMatrix"), s=0)
        #Train using only age and gender
        model_ag[[i]][[j]][[k]]=glmnet(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))),train_cv_input$response_cat, weights=train_cv_input$weight, family='binomial', lambda = 0)
        pred_train_ag[[i]][[j]][[k]]=predict(model_ag[[i]][[j]][[k]], newx=as(t(rbind(as.numeric(train_cv_input$age_imputed), as.numeric(train_cv_input$gender_imputed=='Male'))),"dgCMatrix"), s=0)
        pred_val_ag[[i]][[j]][[k]]=predict(model_ag[[i]][[j]][[k]], newx=as(t(rbind(as.numeric(val_cv_input$age_imputed), as.numeric(val_cv_input$gender_imputed=='Male'))),"dgCMatrix"), s=0)
        pred_test_ag[[i]][[j]][[k]]=predict(model_ag[[i]][[j]][[k]], newx=as(t(rbind(as.numeric(test_input$age_imputed), as.numeric(test_input$gender_imputed=='Male'))),"dgCMatrix"), s=0)
      }
      
      #Compute accuracy/AUC in CV folds
      accuracy_train_cv[k]=length(which((pred_train[[i]][[j]][[k]]>0.5)==train_cv_input$response_cat))/ncol(train_cv_input)
      accuracy_val_cv[k]=length(which((pred_val[[i]][[j]][[k]]>0.5)==val_cv_input$response_cat))/ncol(val_cv_input)
      accuracy_test_cv[k]=length(which((pred_test[[i]][[j]][[k]]>0.5)==test_input$response_cat))/ncol(test_input)
      accuracy_train_ag_cv[k]=length(which((pred_train_ag[[i]][[j]][[k]]>0.5)==train_cv_input$response_cat))/ncol(train_cv_input)
      accuracy_val_ag_cv[k]=length(which((pred_val_ag[[i]][[j]][[k]]>0.5)==val_cv_input$response_cat))/ncol(val_cv_input)
      accuracy_test_ag_cv[k]=length(which((pred_test_ag[[i]][[j]][[k]]>0.5)==test_input$response_cat))/ncol(test_input)
      auc_train_cv[k]=auc(train_cv_input$response_cat,as.vector(pred_train[[i]][[j]][[k]]))
      auc_val_cv[k]=if(length(unique(val_cv_input$response_cat))==2) auc(val_cv_input$response_cat,as.vector(pred_val[[i]][[j]][[k]])) else NA
      auc_test_cv[k]=if(length(unique(test_input$response_cat))==2) auc(test_input$response_cat,as.vector(pred_test[[i]][[j]][[k]])) else NA
      auc_train_ag_cv[k]=auc(train_cv_input$response_cat,as.vector(pred_train_ag[[i]][[j]][[k]]))
      auc_val_ag_cv[k]=if(length(unique(val_cv_input$response_cat))==2) auc(val_cv_input$response_cat,as.vector(pred_val_ag[[i]][[j]][[k]])) else NA
      auc_test_ag_cv[k]=if(length(unique(test_input$response_cat))==2) auc(test_input$response_cat,as.vector(pred_test_ag[[i]][[j]][[k]])) else NA
    }
    #Average accuracy/AUC over CV folds
    accuracy_train[[i]][j]=mean(accuracy_train_cv)
    accuracy_train_std[[i]][j]=sd(accuracy_train_cv)
    accuracy_val[[i]][j]=mean(accuracy_val_cv)
    accuracy_val_std[[i]][j]=sd(accuracy_val_cv)
    accuracy_test[[i]][j]=mean(accuracy_test_cv)
    accuracy_test_std[[i]][j]=sd(accuracy_test_cv)
    accuracy_train_ag[[i]][j]=mean(accuracy_train_ag_cv)
    accuracy_train_ag_std[[i]][j]=sd(accuracy_train_ag_cv)
    accuracy_val_ag[[i]][j]=mean(accuracy_val_ag_cv)
    accuracy_val_ag_std[[i]][j]=sd(accuracy_val_ag_cv)
    accuracy_test_ag[[i]][j]=mean(accuracy_test_ag_cv)
    accuracy_test_ag_std[[i]][j]=sd(accuracy_test_ag_cv)
    auc_train[[i]][j]=mean(auc_train_cv)
    auc_train_std[[i]][j]=sd(auc_train_cv)
    auc_val[[i]][j]=mean(auc_val_cv)
    auc_val_std[[i]][j]=sd(auc_val_cv)
    auc_test[[i]][j]=mean(auc_test_cv)
    auc_test_std[[i]][j]=sd(auc_test_cv)
    auc_train_ag[[i]][j]=mean(auc_train_ag_cv)
    auc_train_ag_std[[i]][j]=sd(auc_train_ag_cv)
    auc_val_ag[[i]][j]=mean(auc_val_ag_cv)
    auc_val_ag_std[[i]][j]=sd(auc_val_ag_cv)
    auc_test_ag[[i]][j]=mean(auc_test_ag_cv)
    auc_test_ag_std[[i]][j]=sd(auc_test_ag_cv)
  }
}

#Plot Results
for (i in 1:length(exp_FC)) {
    df=rbind('AUC.GE'=auc_test[[i]],'Std.GE'=auc_test_std[[i]])
    df=data.frame(t(df))
    df$Vaccine=rownames(df)
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
  ggsave(paste('LOVO_wCV_Test_Accuracies_Day',tp_int[i+1], '_nfeat', n_feat_select, '_', 'lasso', '.pdf',sep=''), width = 7, height = 8)
}


# Get coefficients of lasso model
# coefs <- lapply(1:length(model), function(x)
#                 lapply(1:length(model[[x]]), function(y)
#                   sapply(1:length(model[[x]][[y]]), function(z)
#                     coef(model[[x]][[y]][[z]], best_lam_list[[x]][[y]][[z]]))))
# 
# # Get mean and SD of coeffs across all 10 CV models 
# for(day in 1:length(coefs)) {
#   dfMeanMaster <- data.frame(matrix(vector(), n_feat_select+3, 1))
#   dfMeanMaster <- dfMeanMaster %>% mutate(gene=rownames(coefs[[day]][[1]][[1]]))
#   dfSDMaster <- data.frame(matrix(vector(), n_feat_select+3, 1))
#   dfSDMaster <- dfSDMaster %>% mutate(gene=rownames(coefs[[day]][[1]][[1]]))
#   for(pt in 1:length(coefs[[day]])) {
#     df <- data.frame(matrix(vector(), n_feat_select+3, 1))
#     df <- df %>% mutate(gene=rownames(coefs[[day]][[pt]][[1]]))
#     for(cv in 1:length(coefs[[day]][[pt]])) {
#       coefValues <- as.data.frame(as.matrix(coefs[[day]][[pt]][[cv]]))
#       coefValues <- coefValues %>% mutate(gene=rownames(coefValues))
#       df <- full_join(df, coefValues, by="gene")
#     }
#     rownames(df) <- df$gene
#     df <- subset(df, select=-gene)
#     #Get mean across all 10 CVs
#     dfMean <- rowMeans(df, na.rm = T)
#     dfMean <- as.data.frame(as.matrix(dfMean))
#     dfMean <- dfMean %>% mutate(gene=rownames(dfMean))
#     #Get standard deviation across all 10 CVs
#     dfSD <- apply(df,1, sd, na.rm = TRUE)
#     dfSD <- as.data.frame(as.matrix(dfSD))
#     colnames(dfSD) = paste("sd", pt_uni_tp[[day]][[pt]], sep='_')
#     dfSD <- dfSD %>% mutate(gene=rownames(dfSD))
#     
#     dfMeanMaster <- full_join(dfMeanMaster, dfMean, by="gene", all=T)  #inner_join or full_join?
#     dfSDMaster <- full_join(dfSDMaster, dfSD, by="gene", all=T)
#   }
#   dfMeanMaster <- dfMeanMaster[,-1]
#   dfSDMaster <- dfSDMaster[,-1]
#   colnames(dfMeanMaster) <- c("gene", pt_uni_tp[[day]])
#   dfMeanMaster <- dfMeanMaster[order(-dfMeanMaster[,2]),]
#   fwrite(dfMeanMaster, file = paste("coeffs_", n_feat_select, 'Genes_' , tp_int[day+1],'.csv',sep=''), row.names = F)
#   fwrite(dfSDMaster, file = paste("coeffs_SD_", n_feat_select, 'Genes_' , tp_int[day+1],'.csv',sep=''), row.names = F)
# }
# 
# 
