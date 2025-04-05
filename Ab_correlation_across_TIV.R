# Perform correlation of individual genes with Ab response

library(RColorBrewer)
library(ggthemes)
library(stats)
library(tidyverse)
library(dplyr)
library(psych)
library(metap)
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

rm(list=ls())
setwd("/Users/audrey/Desktop/hipc/")
#Load master color set
load("adj_path_vt_colors.rda")
#Load GE
load("./2020_01_22_IS2_withResponse_noNormCorr.rda")
#Rename to common variable
eset=IS2_withResponse_noNormCorr
rm(IS2_withResponse_noNormCorr)
#Add subject colnames
colnames(eset)=colnames(exprs(eset))
#Cutoff/input options
p_cutoff=0.05
fdr_cutoff=0.05
NES_cutoff=2
corrMethod <- "pearson"

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
eset=eset[,-which(eset$MFC == Inf)]  #exclude study1289 which has MFC = Inf

#Find samples from timepoints of interest
tp_int=c(0,1,3,7)
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
cv <- apply(exprs(eset), 1, function(x) (sd(x)/mean(x)))  #find genes with largest coeff. of variance
top_cv <- order(-cv)[1:1000]

#Create combined SDY/pathogen/vaccine type column
eset$SDY_pathogen_type=paste(eset$study,eset$pathogen_type)
#Create unique list of studies
matrix_uni=unique(eset$matrix)

#Create master color list for pathogen/vaccine type combination
cols=colorRampPalette(brewer.pal(length(unique(eset$pathogen_type)), "Set3"))
mycolors=cols(length(unique(eset$pathogen_type)))
names(mycolors)=unique(eset$pathogen_type)

#Collapse to BTM level
#Convert genes to entrez ID
entrez_gene_lookup=select(org.Hs.eg.db, rownames(eset), c("ENTREZID"), "ALIAS")
#Remove genes with multiple entrez IDs
entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,2]) | duplicated(entrez_gene_lookup[,2], fromLast=TRUE)),]
#Remove entrez IDs with multiple genes
entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,1]) | duplicated(entrez_gene_lookup[,1], fromLast=TRUE)),]
entrez_list=as.character(entrez_gene_lookup[match(rownames(eset),entrez_gene_lookup[,1]),2])  #get entrezID of genes that are in eset
#Load BTMs
load('./BTM_for_GSEA_20131008_geneID.RData')
#Remove BTMs which have no matching genes in dataset
BTM_list=BTM_list[!sapply(BTM_list, function(x) is_empty(intersect(x,entrez_list)))]
#Collapse - match entrezID in BTM_list to entrezID of genes (arrange according to rownames of eset) - find arithmetic mean for each sample
exp_BTM=do.call(rbind, lapply(BTM_list, function(x) colMeans(exprs(eset[na.omit(match(x,entrez_list)),]),na.rm=TRUE)))
#Create BTM eset (maybe not proper approach but it works)
eset_BTM=eset[1:nrow(exp_BTM),]
rownames(eset_BTM)=rownames(exp_BTM)
exprs(eset_BTM)=exp_BTM

#Compute D0 normalized FC - for individual genes - exp_FC has 3 large matrix with D1, D3, D7 timepoints
ind_D0=which(0==eset$study_time_collected)
common=lapply(2:length(ind),function(x) intersect(eset$participant_id[ind[[x]]],eset$participant_id[ind_D0]))  #get participant id with timepoint and D0 data
ia=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind[[x]]]))) #included participant IDs match with participant IDs at diff timepoints in eset
ib=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind_D0])))  #participant IDs match with D0 to get indices in ind_D0
exp_FC=lapply(2:length(ind),function(x) eset[,ind[[x]][ia[[x-1]]]])  #get expression data of the 3 timepoints
exp_FC=lapply(2:length(ind),function(x) {exprs(exp_FC[[x-1]])=exprs(exp_FC[[x-1]])-exprs(eset[,ind_D0[ib[[x-1]]]]); exp_FC[[x-1]]}) 
#Compute D0 normalized FC for BTMs
exp_BTM_FC=lapply(2:length(ind),function(x) eset_BTM[,ind[[x]][ia[[x-1]]]])
exp_BTM_FC=lapply(2:length(ind),function(x)
{exprs(exp_BTM_FC[[x-1]])=exprs(exp_BTM_FC[[x-1]])-exprs(eset_BTM[,ind_D0[ib[[x-1]]]]); exp_BTM_FC[[x-1]]})

#For each study, average expression across all subjects
matrix_uni_tp=lapply(exp_FC,function(x) x$matrix[!duplicated(x$matrix)])  #study matrix for each timepoint
#Store study metadata
matrix_uni_tp_metaData=lapply(exp_FC,function(x) pData(x)[!duplicated(x$matrix),] %>%
                                rownames_to_column() %>%
                                dplyr::select(matrix, vaccine, vaccine_type, pathogen, adjuvant, pathogen_type) %>%
                                column_to_rownames(var = "matrix"))

#Study size
#exp_FC_n is the study size, matrix_uni_tp is the corresponding study, matrix_ind is the indices of the corresponding study
matrix_ind=lapply(1:length(exp_FC),
                  function(x) lapply(1:length(matrix_uni_tp[[x]]), function(y) which(matrix_uni_tp[[x]][[y]]==exp_FC[[x]]$matrix)))
exp_FC_n=lapply(1:length(exp_FC),
                function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) length(matrix_ind[[x]][[y]])))

#Average FC
exp_FC_mean=lapply(1:length(exp_FC),
                   function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) rowMeans(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]))))
exp_BTM_FC_mean=lapply(1:length(exp_BTM_FC),
                       function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) rowMeans(exprs(exp_BTM_FC[[x]][,matrix_ind[[x]][[y]]]))))


#cor.test for individual genes - each timepoint, for each study, for each study dataframe row - AL
exp_FC_corr_less=
  lapply(1:length(exp_FC),  #for each timepoint
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  #for each study
                            function(y) apply(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  #for each gene
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "less")$p.value)))
exp_FC_corr_greater=
  lapply(1:length(exp_FC),  
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  
                            function(y) apply(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "greater")$p.value)))

#cor.test for BTM
exp_BTM_FC_corr_less=
  lapply(1:length(exp_BTM_FC),  #for each timepoint
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  #for each study
                            function(y) apply(exprs(exp_BTM_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  #for each BTM
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "less")$p.value)))
exp_BTM_FC_corr_greater=
  lapply(1:length(exp_BTM_FC),  
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  
                            function(y) apply(exprs(exp_BTM_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "greater")$p.value)))

##Get estimate (correlation coeff r) - two-sided cor.test
exp_FC_corr_estimate=
  lapply(1:length(exp_FC),  #for each timepoint
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  #for each study
                            function(y) apply(exprs(exp_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  #for each gene
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "two.sided")$estimate)))
exp_BTM_FC_corr_estimate=
  lapply(1:length(exp_BTM_FC),  
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  
                            function(y) apply(exprs(exp_BTM_FC[[x]][,matrix_ind[[x]][[y]]]), 1,  
                                              function(z) cor.test(z, pData(exp_FC[[x]])[matrix_ind[[x]][[y]],"MFC"], method = corrMethod, alternative = "two.sided")$estimate)))


#Replace NA p values with 0.9999 
#cor.test gives NA, t.test gives NaN
exp_FC_corr_less=lapply(exp_FC_corr_less, function(x) replace(x,is.na(x),0.9999))
exp_FC_corr_greater=lapply(exp_FC_corr_greater, function(x) replace(x,is.na(x),0.9999))
exp_BTM_FC_corr_greater = lapply(exp_BTM_FC_corr_greater, function(x) replace(x,is.na(x),0.9999))
exp_BTM_FC_corr_less = lapply(exp_BTM_FC_corr_less, function(x) replace(x,is.na(x),0.9999))

#Replace R values = 1 with 0.999. NA r values are omitted later when finding mean
exp_FC_corr_estimate=lapply(exp_FC_corr_estimate, function(x) replace(x,x==1,0.9999))
exp_BTM_FC_corr_estimate=lapply(exp_BTM_FC_corr_estimate, function(x) replace(x,x==1,0.9999))

#Replace 1 p values with 0.9999 (for compatibility with sumz)
exp_FC_corr_less=lapply(exp_FC_corr_less, function(x) replace(x,x>0.9999,0.9999))
exp_FC_corr_greater=lapply(exp_FC_corr_greater, function(x) replace(x,x>0.9999,0.9999))
exp_BTM_FC_corr_greater=lapply(exp_BTM_FC_corr_greater, function(x) replace(x,x>0.9999,0.9999))
exp_BTM_FC_corr_less=lapply(exp_BTM_FC_corr_less, function(x) replace(x,x>0.9999,0.9999))

#Find unique pathogen/vaccine types by timepoint
pathogen_type_uni_tp=lapply(matrix_uni_tp_metaData,function(x) x$pathogen_type[!duplicated(x$pathogen_type)])
pathogen_type_metaData=lapply(1:length(matrix_uni_tp_metaData),
                              function(x) dplyr::select(
                                matrix_uni_tp_metaData[[x]][!duplicated(matrix_uni_tp_metaData[[x]]$pathogen_type),],
                                vaccine_type, pathogen, adjuvant))


# Run Stouffer's method for the pvals for each pathogen type - corrected pvals are stored in a list of ea timepoint having a df with row=genes, col=pathogentype
pathogen_type_ind=vector("list",length(exp_FC))
pathogen_type_FC_mean=vector("list",length(exp_FC))
pathogen_type_p=vector("list",length(exp_FC)) #to store integrated pval after Stouffer's method
pathogen_type_q=vector("list",length(exp_FC))
pathogen_type_fisherz=vector("list",length(exp_FC)) #store fisher Z 
pathogen_type_fisherz_Rvalues=vector("list",length(exp_FC))
pathogen_type_fisherz_mean=vector("list",length(exp_FC))
pathogen_type_up_down_logic=vector("list",length(exp_FC))
pathogen_type_DEG_up=vector("list",length(exp_FC))
pathogen_type_DEG_down=vector("list",length(exp_FC))
pathogen_type_DEG_total_num=vector("list",length(exp_FC))
pathogen_type_DEG_up_overlap=vector("list",length(exp_FC))
pathogen_type_DEG_down_overlap=vector("list",length(exp_FC))

BTM_pathogen_type_FC_mean=vector("list",length(exp_FC))
BTM_pathogen_type_p=vector("list",length(exp_FC)) #to store integrated pval after Stouffer's method
BTM_pathogen_type_q=vector("list",length(exp_FC))
BTM_pathogen_type_fisherz=vector("list",length(exp_FC)) #store fisher Z 
BTM_pathogen_type_fisherz_Rvalues=vector("list",length(exp_FC))
BTM_pathogen_type_fisherz_mean=vector("list",length(exp_FC))
BTM_pathogen_type_up_down_logic=vector("list",length(exp_FC))
BTM_pathogen_type_DEG_up=vector("list",length(exp_FC))
BTM_pathogen_type_DEG_down=vector("list",length(exp_FC))
BTM_pathogen_type_DEG_total_num=vector("list",length(exp_FC))
BTM_pathogen_type_DEG_up_overlap=vector("list",length(exp_FC))
BTM_pathogen_type_DEG_down_overlap=vector("list",length(exp_FC))

comb=vector("list",length(exp_FC))
for (i in 1:length(exp_FC)) {
  #For each pathogen/type combo
  pathogen_type_FC_mean[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(pathogen_type_uni_tp[[i]]))
  pathogen_type_p[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(pathogen_type_uni_tp[[i]]))
  pathogen_type_q[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(pathogen_type_uni_tp[[i]]))
  pathogen_type_up_down_logic[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(pathogen_type_uni_tp[[i]]))
  pathogen_type_ind[[i]]=vector("list",length(pathogen_type_uni_tp[[i]]))
  pathogen_type_DEG_up[[i]]=vector("list",length(pathogen_type_uni_tp[[i]]))
  pathogen_type_DEG_down[[i]]=vector("list",length(pathogen_type_uni_tp[[i]]))
  pathogen_type_DEG_total_num[[i]]=matrix(NA, nrow=1, ncol=length(pathogen_type_uni_tp[[i]]))
  
  BTM_pathogen_type_FC_mean[[i]]=matrix(NA, nrow=nrow(exp_BTM), ncol=length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_p[[i]]=matrix(NA, nrow=nrow(exp_BTM), ncol=length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_q[[i]]=matrix(NA, nrow=nrow(exp_BTM), ncol=length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_up_down_logic[[i]]=matrix(NA, nrow=nrow(exp_BTM), ncol=length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_DEG_up[[i]]=vector("list",length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_DEG_down[[i]]=vector("list",length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_DEG_total_num[[i]]=matrix(NA, nrow=1, ncol=length(pathogen_type_uni_tp[[i]]))
  
  pathogen_type_fisherz_mean[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(pathogen_type_uni_tp[[i]]))
  BTM_pathogen_type_fisherz_mean[[i]]=matrix(NA, nrow=nrow(exp_BTM), ncol=length(pathogen_type_uni_tp[[i]]))
  
  for (j in 1:length(pathogen_type_uni_tp[[i]])) {
    pathogen_type_ind[[i]][[j]]=which(pathogen_type_uni_tp[[i]][[j]]==matrix_uni_tp_metaData[[i]]$pathogen_type)  #get indices of individual pathogen type in the meta table
    #If there is more than 1 study, integrate p values by Stouffer's method
    if (length(pathogen_type_ind[[i]][[j]])>1) {  #if there are > 1 study for that pathogen type
      #Average FCs (weighted by n)
      pathogen_type_FC_mean[[i]][,j]=rowWeightedMeans(exp_FC_mean[[i]][,pathogen_type_ind[[i]][[j]]], w=exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]])
      BTM_pathogen_type_FC_mean[[i]][,j]=rowWeightedMeans(exp_BTM_FC_mean[[i]][,pathogen_type_ind[[i]][[j]]], w=exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]])
      #Directional Stouffer's method using sumz function (weighted by sqrt(n))
      pathogen_type_p_less=vapply(1:nrow(eset), FUN.VALUE=1, function(x)
        sumz(exp_FC_corr_less[[i]][x,pathogen_type_ind[[i]][[j]]],sqrt(exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]]))$p)
      pathogen_type_p_greater=vapply(1:nrow(eset), FUN.VALUE=1, function(x)
        sumz(exp_FC_corr_greater[[i]][x,pathogen_type_ind[[i]][[j]]],sqrt(exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]]))$p)
      pathogen_type_up_down_logic[[i]][,j]=pathogen_type_p_greater < pathogen_type_p_less #logical variable, true if upregulation p value is smallest
      pathogen_type_p[[i]][,j]=2*pmin(pathogen_type_p_less,pathogen_type_p_greater)
      
      #for BTM
      #Average FCs (weighted by n)
      BTM_pathogen_type_FC_mean[[i]][,j]=rowWeightedMeans(exp_BTM_FC_mean[[i]][,pathogen_type_ind[[i]][[j]]], w=exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]])
      #Directional Stouffer's method using sumz function (weighted by sqrt(n))
      BTM_pathogen_type_p_less=vapply(1:nrow(exp_BTM), FUN.VALUE=1, function(x)
        sumz(exp_BTM_FC_corr_less[[i]][x,pathogen_type_ind[[i]][[j]]],sqrt(exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]]))$p)
      BTM_pathogen_type_p_greater=vapply(1:nrow(exp_BTM), FUN.VALUE=1, function(x)
        sumz(exp_BTM_FC_corr_greater[[i]][x,pathogen_type_ind[[i]][[j]]],sqrt(exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]]))$p)
      BTM_pathogen_type_up_down_logic[[i]][,j]= BTM_pathogen_type_p_greater < BTM_pathogen_type_p_less #logical variable, true if upregulation p value is smallest
      BTM_pathogen_type_p[[i]][,j]=2*pmin(BTM_pathogen_type_p_less, BTM_pathogen_type_p_greater) 
      
      #Integrate pearson correlation across studies of same pathogen type - for each gene
      pathogen_type_fisherz[[i]] = apply(exp_FC_corr_estimate[[i]], 2, function(x) FisherZ(x))
      BTM_pathogen_type_fisherz[[i]] = apply(exp_BTM_FC_corr_estimate[[i]], 2, function(x) FisherZ(x))
      #Get weighted mean of Fisher Z for each pathogen type 
      pathogen_type_fisherz_mean[[i]][,j]=rowWeightedMeans(pathogen_type_fisherz[[i]][,pathogen_type_ind[[i]][[j]]], w=exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]], na.rm = T)
      BTM_pathogen_type_fisherz_mean[[i]][,j]=rowWeightedMeans(BTM_pathogen_type_fisherz[[i]][,pathogen_type_ind[[i]][[j]]], w=exp_FC_n[[i]][pathogen_type_ind[[i]][[j]]], na.rm = T)
      
      #Otherwise store single p values
    } else {
      pathogen_type_FC_mean[[i]][,j]=exp_FC_mean[[i]][,pathogen_type_ind[[i]][[j]]]
      pathogen_type_p[[i]][,j]=2*pmin(exp_FC_corr_less[[i]][,pathogen_type_ind[[i]][[j]]],exp_FC_corr_greater[[i]][,pathogen_type_ind[[i]][[j]]])
      pathogen_type_up_down_logic[[i]][,j]=exp_FC_corr_greater[[i]][,pathogen_type_ind[[i]][[j]]] < exp_FC_corr_less[[i]][,pathogen_type_ind[[i]][[j]]]
      
      #for BTM
      BTM_pathogen_type_FC_mean[[i]][,j]=exp_BTM_FC_mean[[i]][,pathogen_type_ind[[i]][[j]]]
      BTM_pathogen_type_p[[i]][,j]=2*pmin(exp_BTM_FC_corr_less[[i]][,pathogen_type_ind[[i]][[j]]],exp_BTM_FC_corr_greater[[i]][,pathogen_type_ind[[i]][[j]]])
      BTM_pathogen_type_up_down_logic[[i]][,j]=exp_BTM_FC_corr_greater[[i]][,pathogen_type_ind[[i]][[j]]] < exp_BTM_FC_corr_less[[i]][,pathogen_type_ind[[i]][[j]]]
      
      #for R correlation coeff
      pathogen_type_fisherz_mean[[i]][,j]=pathogen_type_fisherz[[i]][,pathogen_type_ind[[i]][[j]]]
      BTM_pathogen_type_fisherz_mean[[i]][,j]=BTM_pathogen_type_fisherz[[i]][,pathogen_type_ind[[i]][[j]]]
    }
    #Correct p values >1
    pathogen_type_p[[i]][(pathogen_type_p[[i]][,j]>1),j]=1
    BTM_pathogen_type_p[[i]][(BTM_pathogen_type_p[[i]][,j]>1),j]=1
    
    #Convert p values to q values (FDR correction) - for each pathogen type
    pathogen_type_q[[i]][,j]=p.adjust(pathogen_type_p[[i]][,j],method="BH") #BH method
    BTM_pathogen_type_q[[i]][,j]=p.adjust(BTM_pathogen_type_p[[i]][,j],method="BH") #BH method
    
    #Convert average Fisher Z for each pathogen type to pearson R coeff
    pathogen_type_fisherz_Rvalues[[i]] = apply(pathogen_type_fisherz_mean[[i]], 2, function(x) FisherZInv(x))
    BTM_pathogen_type_fisherz_Rvalues[[i]] = apply(BTM_pathogen_type_fisherz_mean[[i]], 2, function(x) FisherZInv(x))
    
    #Find DEGs
    ind_up=which(pathogen_type_up_down_logic[[i]][,j]==TRUE)
    ind_down=which(pathogen_type_up_down_logic[[i]][,j]==FALSE)
    ind_sig=which(pathogen_type_q[[i]][,j] < fdr_cutoff)
    pathogen_type_DEG_up[[i]][[j]]=intersect(ind_up,ind_sig)
    pathogen_type_DEG_down[[i]][[j]]=intersect(ind_down,ind_sig)
    pathogen_type_DEG_total_num[[i]][j]=length(ind_sig)
    
    #Find significant BTMs
    BTM_ind_up=which(BTM_pathogen_type_up_down_logic[[i]][,j]==TRUE)
    BTM_ind_down=which(BTM_pathogen_type_up_down_logic[[i]][,j]==FALSE)
    BTM_ind_sig=which(BTM_pathogen_type_q[[i]][,j] < fdr_cutoff)
    BTM_pathogen_type_DEG_up[[i]][[j]]=intersect(BTM_ind_up,BTM_ind_sig)
    BTM_pathogen_type_DEG_down[[i]][[j]]=intersect(BTM_ind_down,BTM_ind_sig)
    BTM_pathogen_type_DEG_total_num[[i]][j]=length(BTM_ind_sig)
  }
  colnames(pathogen_type_FC_mean[[i]])=pathogen_type_uni_tp[[i]]
  rownames(pathogen_type_FC_mean[[i]])=rownames(fData(eset))
  colnames(pathogen_type_p[[i]])=pathogen_type_uni_tp[[i]]
  colnames(pathogen_type_q[[i]])=pathogen_type_uni_tp[[i]]
  rownames(pathogen_type_metaData[[i]])=pathogen_type_uni_tp[[i]]
  
  colnames(BTM_pathogen_type_FC_mean[[i]])=pathogen_type_uni_tp[[i]]
  rownames(BTM_pathogen_type_FC_mean[[i]])=rownames(exp_BTM)
  colnames(BTM_pathogen_type_p[[i]])=pathogen_type_uni_tp[[i]]
  colnames(BTM_pathogen_type_q[[i]])=pathogen_type_uni_tp[[i]]
  
  colnames(pathogen_type_fisherz_Rvalues[[i]])=pathogen_type_uni_tp[[i]]
  colnames(BTM_pathogen_type_fisherz_Rvalues[[i]])=pathogen_type_uni_tp[[i]]
  rownames(pathogen_type_fisherz_Rvalues[[i]])=rownames(fData(eset))
  rownames(BTM_pathogen_type_fisherz_Rvalues[[i]])=rownames(exp_BTM)
}

#Look at within TIV studies 
vaccine_of_interest <- "Influenza  ( Inactivated )"
tiv_ind=vector("list",length(exp_FC))
exp_FC_corr_p = vector("list",length(exp_FC))
exp_FC_corr_q = vector("list",length(exp_FC))
exp_FC_r = vector("list",length(exp_FC))
BTM_exp_FC_r = vector("list",length(exp_FC))
exp_FC_logic = vector("list",length(exp_FC))
for (i in 1:length(exp_FC)) {
  tiv_ind[[i]] = which(matrix_uni_tp_metaData[[i]]$pathogen_type == vaccine_of_interest)
  study_matrix <- rownames(matrix_uni_tp_metaData[[i]])[tiv_ind[[i]]]
  exp_FC_corr_p[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(tiv_ind[[i]]))
  exp_FC_corr_q[[i]]=matrix(NA, nrow=nrow(eset), ncol=length(tiv_ind[[i]]))
  exp_FC_logic[[i]] = matrix(NA, nrow=nrow(eset), ncol=length(tiv_ind[[i]]))
  for (j in 1:length(tiv_ind[[i]])) {
    #Convert p values to q values (FDR correction)
    exp_FC_corr_p[[i]][,j] = 2*pmin(exp_FC_corr_less[[i]][,tiv_ind[[i]][[j]]],exp_FC_corr_greater[[i]][,tiv_ind[[i]][[j]]])
    exp_FC_logic[[i]][,j] = exp_FC_corr_greater[[i]][,tiv_ind[[i]][[j]]] < exp_FC_corr_less[[i]][,tiv_ind[[i]][[j]]]
    exp_FC_corr_q[[i]][,j] =  p.adjust(exp_FC_corr_p[[i]][,j],method="BH") #BH method
    exp_FC_r[[i]] = exp_FC_corr_estimate[[i]][,tiv_ind[[i]]]
    BTM_exp_FC_r[[i]] = exp_BTM_FC_corr_estimate[[i]][,tiv_ind[[i]]]
    new_study_matrix <- sapply(study_matrix, function(x) str_remove(x, "WHOLEBLOOD_"))  #get only SDY for names
    new_study_matrix <- sapply(new_study_matrix, function(x) str_remove(x, "PBMC_"))
    new_study_matrix <- sapply(new_study_matrix, function(x) str_remove(x, "_GEO"))
    new_study_matrix <- sapply(new_study_matrix, function(x) str_remove(x, "_IMMPORT"))
    colnames(exp_FC_r[[i]]) <- new_study_matrix
    colnames(BTM_exp_FC_r[[i]]) <- new_study_matrix
  }
}
ind_sig_genes_tiv <- lapply(exp_FC_corr_q, function(y) apply(y, 2, function(x) which(x < 0.1)))

#for BTM
colorLS=colorRampPalette(colors = c("blue", "white", "red"))(n = 100)
BTM_tiv_pairwise_corr <- lapply(BTM_exp_FC_r, cor, method="spearman", use="complete.obs")
corrplot(BTM_tiv_pairwise_corr[[3]], method="circle", col=colorLS, tl.cex = 0.7, tl.srt = 45, tl.col='black', win.asp = 1)