## ---------------------------
##
## Script name: Foundin_main.R
##
## Purpose of script: Includes: 1) a function to generate correlation of expression of a subset of genes to the differentiation efficiency as estimated by percentage of cells of certain 
## subtype 
##				2) Set of commands to run to get predictive of differentiation genes based on their expression table and diff.efficiency annotation(binary; should be factor, e.g. "Y"/"N" (high diff/efficiency?-Yes/No)
##				3) Generation of diagnostic plots for predictive model
## Author: Natalia Savytska
##
## Date Created: 2021-03-24, last edited 2021-05-27
##
## Copyright (c) Natalia Savytska, 2021
## Email: savytskanat@gmail.com



library(glmnet)
library(caret)
library(data.table)
library(ggplot2)
library(ROCR)



####################################################################################################################
####################################################################################################################


## Function to make a correlation analysis


# As input corr_efficiency() takes:
#                       -ExprCnt - data.frame, where columns are samples (with respective sample IDs as colnames) and rows are normalized counts per gene (geneIDs as rownames)
#                       -CellCnt - data.frame, where columns are Cell types and rows are counts for each cell type per Sample (rownames should be sample IDs; sample ID conventions should be same as in ExprCnt colnames)
#                       -gene_id - character vector, that contains IDs of genes one is interested to correlate. Should be same ID convention as in rownames of ExprCnt
#                       -cell_id - character string, that contains column name for cell type one is interested to correlate Expression of genes to; should be one of colnames in CellCnt
#                       -corr_method - character string, indicating method for correlation test: "kendall" (default), "spearman" or "pearson"
# As output returns data.frame with correlation coefficient, pval and adjusted pval (FDR) for each gene of interest in gene_id vector
  


corr_efficiency<-function(ExprCnt,CellCnt,gene_id,cell_id="Dopaminergic Neurons",corr_method="kendall"){

  ExprCnt<-ExprCnt[,colnames(ExprCnt) %in% rownames(CellCnt)]
  #order(match(names(dds_zhang_nc_naive_zs),zhang_meta2$a))
  ExprCnt<-ExprCnt[,order(match(colnames(ExprCnt),rownames(CellCnt)))]
  if (identical(rownames(CellCnt),colnames(ExprCnt))==FALSE){
    print("Something is off with sample IDs in either Expression Table or Cell Count Table.\n
          Please, check the ID conventions in both. They should have common IDs")
  } else if (identical(rownames(CellCnt),colnames(ExprCnt))==TRUE){
    ExprCnt_T<-as.data.frame(t(ExprCnt))
    
    df_coeff_g<-data.frame(Genes=gene_id,corrval=0,pval=0)
    for (i in gene_id){
      df_coeff_g[df_coeff_g$Genes==i,"corrval"]<-cor.test(ExprCnt_T[,i],CellCnt[,cell_id],method=corr_method,na.rm=TRUE)$estimate
      df_coeff_g[df_coeff_g$Genes==i,"pval"]<-cor.test(ExprCnt_T[,i],CellCnt[,cell_id],method=corr_method,na.rm=TRUE)$p.value
    }
    df_coeff_g$FDR<-p.adjust(df_coeff_g$pval,method="fdr")
    return(df_coeff_g)
  } else {
    print("Something is off with sample IDs or gene IDs or format of the tables")
  }
}
  

# Plot histos for obrained unadjusted pvals and coefficients
    
ggplot(corr_df, aes(x = pval)) +
      geom_histogram(position = "identity", alpha = .8, bins=40)
    
ggplot(corr_df, aes(x = corrval)) +
      geom_histogram(position = "identity", alpha = .8, bins=40)




####################################################################################################################
####################################################################################################################





# LOGISTIC REGRESSION, tuning+final model


# GLM


## Load t_counts
transposed_normalizedCount_table <- read.delim("transposed_normalizedCount_table.txt", row.names=1)

# Relevel
transposed_normalizedCount_table[,16397]<-as.factor(transposed_normalizedCount_table[,16397])
transposed_normalizedCount_table[,16397]<-relevel(transposed_normalizedCount_table[,16397],"Y")


# As input, one should provide a transposed count table (tCntTab), in which Samples/observations are in rows, and Genes/predictors are in columns.
# This tCntTable should also contain a (factor/character) column with the feature to be predicted, in this script - DiffEfficiencyColID (in our tables values for this are "Y" and "N"; can be "High" and "Low" or anything else).



# Training model parameters.
# Method: "cv" for N-fold crossvalidation; number for how many folds (N) data will be split into; seeds - seeds
# Generate seed numbers for random sampling in crossvalidations
set.seed(25)
seeds_glm <- vector(mode = "list", length = 6)
for(i in 1:5) seeds_glm[[i]] <- sample.int(1000, 3)

## For the last model:
set.seed(25)
seeds_glm[[6]] <- sample.int(1000, 1)

trctrl <- trainControl(method = "cv",number=5, seeds=seeds_glm)


# First run wih default setting, to detect best lambda
set.seed(25)
enetFit_dopa_el <- train(DopaNeu_High~., data = transposed_normalizedCount_table, 
                         method = "glmnet",
                         trControl=trctrl,
                         # alpha and lambda paramters to try
                         tuneGrid = expand.grid(.alpha = 0.5))
# Resulting best lambda: 0.22; 

# alpha lambda  Accuracy     Kappa AccuracySD   KappaSD
# 0.5   0.22 0.7827614 0.1872604 0.03434338 0.1733364


# Run with set lambda

set.seed(25)
enetFit_dopa_el <- train(DopaNeu_High~., data = transposed_normalizedCount_table, 
                         method = "glmnet",
                         trControl=trctrl,
                         # alpha and lambda paramters to try
                         tuneGrid = expand.grid(.alpha = 0.5, .lambda=0.22))



# Assess the trained model
enetFit_dopa_el


####################################################################################################################
####################################################################################################################

# COEFFICIENT PLOT



# Plotting features` (non-zero) coefficients
# Plot Coefficients for the Candidate Predictors 
coeff<-as.data.frame((as.matrix(coef(enetFit_dopa_el$finalModel,enetFit_dopa_el$bestTune$lambda))))
coeff$ID<-rownames(coeff)
colnames(coeff)<-c("Coeff","ID")
coeff<-coeff[coeff$Coeff!=0,]


coeff_cand<-coeff[coeff$Coeff!=0,]
coeff_cand<-coeff_cand[-1,]
coeff_cand<-coeff_cand[order(coeff_cand$Coeff),]
coeff_cand$GeneID<-factor(coeff_cand$GeneID,levels = coeff_cand$GeneID)

# Plot
ggplot(coeff_cand,aes(x=GeneID, y = Coeff))+ 
  geom_bar(stat="identity", position="identity")+ xlab("")+ylab("Coefficients")+
  geom_text(aes(label=c(round(coeff_cand$Coeff, digits4)[1:4],"","","","","","")), position=position_dodge(width=0), vjust=0.5, hjust=1, size = 8, colour="grey31")+
  geom_text(aes(label=c("","","","",round(coeff_cand$Coeff, digits4)[5:10])), position=position_dodge(width=0), vjust=0.5, hjust=0, size = 8, colour="grey31")+
  coord_flip() + theme_minimal() +
  theme(axis.text.x = element_text( hjust = 1, vjust = 0.5, size = 20),
        axis.text.y = element_text( hjust = 1, vjust = 0.5, size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_text(size=25,face="bold",margin = margin(t = 20, r = 20, b = 0, l = 20)))  + scale_y_continuous(limits = c(-0.012, 0.01), breaks =c(-0.01,-0.005,0.0,0.005,0.01))




####################################################################################################################
####################################################################################################################

#CONFUSION MX PLOT



# Plotting confusion mx

## Plot Confusion mx
# Get a mx for plotting
# Get probabilities
class.res_f<-predict(enetFit_dopa_el,transposed_normalizedCount_table[,-16397],type="prob")
conf_mx<-data.frame(True=transposed_normalizedCount_table[,16397],Predicted=class.res_f[,"Y"])
# I want all False hits red, all True hits Blue
# red - #CC0033 ; teal #3399CC
# Optimal cutoff 0.74
conf_mx$Color<-""
conf_mx[(conf_mx$True=="Y" & conf_mx$Predicted>0.74) | (conf_mx$True=="N" & conf_mx$Predicted<=0.74),"Color"]<-"#3399CC"
conf_mx[(conf_mx$True=="Y" & conf_mx$Predicted<=0.74) | (conf_mx$True=="N" & conf_mx$Predicted>0.74),"Color"]<-"#CC0033"


# Plot Confusion mx
ggplot(conf_mx, aes(x=True, y=Predicted, fill=Color))+ 
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=0.75,position="identity")+scale_fill_manual(values=c("#3399CC","#CC0033"))+
  geom_hline(yintercept=0.74, linetype="dashed", color = "grey")+
  scale_x_discrete(labels=c("High","Low"))+
  labs(x="Differentiation Efficiency", y="Predicted Probability")+theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size=20))


####################################################################################################################
####################################################################################################################


# PLOT ROC

# Plot ROC with AUC and other parameters



# Fancy ROC plot
# Plotting ROC and Precision\Recall 
# Get probabilities
class.res_f<-predict(enetFit_dopa_el,transposed_normalizedCount_table[,-16397],type="prob")
pred_d <- prediction(class.res_f[,"Y"], transposed_normalizedCount_table[,16397])
perfm_d <- performance(pred_d, "tpr", "fpr")
auc.perf <- performance(pred_d, measure = "auc")
auc.perf <- unlist(slot(auc.perf, "y.values"))



# Calculate stats for the cutoff 0.74
# Accuracy
pred_vec<-factor(ifelse(class.res_f[,"Y"]>0.74, "Y", "N"),levels=c("Y","N"))
Acc<-round(confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$overall[1],digits=3)
# Accuracy 
# 0.8674699 
Sens<-round(confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[1], digits=3)
# Sensitivity 
# 0.8387097



# FPR
# FPR = FP/(FP+TN) = 1 - Specificity 
# Specificity 
confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[2]
# 0.952381 
FPR<-round((1-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[2]),digits=3)
# 0.04761905

# Precision confusionMatrix(pred_vec,truval)$byClass[5]
Prec<-round(confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[5], digits=3)

# Confidence Interval 95% CI : (0.7752, 0.9319)
# To print upper and lower bound of CI (AccuracyLower, AccuracyUpper) :
ACC_CI<-round(confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$overall[3:4], digits=4)
# AccuracyLower AccuracyUpper 
# 0.7752283     0.9319426


# Plot ROC

plot(perfm_d, main=paste0(("ROC Curve, AUC="),as.character(round(auc.perf,digits=3))), col="red", lwd=3) + abline(a = 0, b = 1, lty=2, lwd=1.5, col="grey") +
  abline(v = FPR, lty=2, lwd=1.5, col="grey31")+abline(h = Sens, lty=2, lwd=1.5, col="grey31")+
  legend(0.6, 0.5,legend=c(paste0("AUC = ",as.character(round(auc.perf,digits=3))), "Optimal Cutoff 0.74",paste0("Accuracy = ",as.character(Acc)), "95% CI : (0.7752, 0.9319)",paste0("TPR = ", Sens), paste0("FPR = ",FPR), paste0("Precision = ", Prec)),
         col = "white", lty=3, bty = "n", cex=1.4)





## Repeated CV (100 times,5-fold). Courtesy: Bansal Vikas, 2021


# Subset to predictive features+ predicted feature
dataX <- transposed_normalizedCount_table[,colnames(transposed_normalizedCount_table) %in% c(coeff$ID[-1],DopaNeu_High)]
auc_vals <- list()
for(i in 1:100){
elasticnet <- cv.glmnet(dataX[,-ncol(dataX)],
                        as.numeric(dataX$DopaNeu_High), family = "binomial", 
                        type.measure = "auc", nfolds=5,  alpha=0.5,lambda = c(0.22,0.24), keep=T)
auc_vals[[i]] <- elasticnet$cvm[2]}summary(unlist(auc_vals))
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.7690  0.8356  0.8565  0.8538  0.8761  0.9124 
sd(unlist(auc_vals))
[1] 0.03134115















####################################################################################################################
####################################################################################################################


## Choosing optimal probability cutoff for predictions

# Get the best cutoff for the model
# Logical vector to numeric with 0`s
num_conv<-function(vec){
  vec<-as.numeric(vec)
  vec<-0
  return(vec)
}

# Create an empty DF with nrow=increments of 0.01 between 0.5 and 0.9; rownames indicate incriment; colnames=c("Sensitivity","Specificity","Precision","Recall","F1","Accuracy")
diffCO_finalModel<-as.data.frame(matrix(,nrow = length(seq(0.5, 0.9, by = 0.01)), ncol = 6))
colnames(diffCO_finalModel)<-c("Sensitivity","Specificity","Precision","Recall","F1","Accuracy")
rownames(diffCO_finalModel)<-paste0("CO_",seq(50, 90, by = 1))


diffCO_finalModel$Sensitivity<-num_conv(diffCO_finalModel$Sensitivity)
diffCO_finalModel$Specificity<-num_conv(diffCO_finalModel$Specificity)
diffCO_finalModel$Precision<-num_conv(diffCO_finalModel$Precision)
diffCO_finalModel$Recall<-num_conv(diffCO_finalModel$Recall)
diffCO_finalModel$F1<-num_conv(diffCO_finalModel$F1)
diffCO_finalModel$Accuracy<-num_conv(diffCO_finalModel$Accuracy)

pred_df<-predict(enetFit_dopa_el,transposed_normalizedCount_table[,-16397], type = "prob" )
# Get the incriment list from 0.5 to 0.9
my_inc<-seq(0.5, 0.9, by = 0.01)
# Fill in the table
for (i in 1:length(seq(50, 90, by = 1))){
  c_inc<-my_inc[i]
  # ifelse(DF_pred[[1]]>0.9, "Y", "N")
  pred_vec<-ifelse(pred_df$Y>c_inc, "Y", "N")
  pred_vec<-as.factor(pred_vec)
  pred_vec<-relevel(pred_vec,"Y")
  
  diffCO_finalModel$Sensitivity[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[1]
  diffCO_finalModel$Specificity[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[2]
  diffCO_finalModel$Precision[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[5]
  diffCO_finalModel$Recall[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[6]
  diffCO_finalModel$F1[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$byClass[7]
  diffCO_finalModel$Accuracy[i]<-confusionMatrix(pred_vec,transposed_normalizedCount_table[,16397])$overall[1]
}
diffCO_finalModel<-diffCO_finalModel[order(diffCO_finalModel$Sensitivity,diffCO_finalModel$Specificity),]

# Get the value for which Sensitivity+Specificity are Max
max(rowSums(diffCO_finalModel[,1:2]))
diffCO_finalModel[rowSums(diffCO_finalModel[,1:2])==max(rowSums(diffCO_finalModel[,1:2])),]
# Best cutoff 0.74
# Sensitivity Specificity Precision    Recall        F1  Accuracy
# CO_74   0.8387097    0.952381 0.9811321 0.8387097 0.9043478 0.8674699

