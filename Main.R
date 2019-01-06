#This script loads all the libraries and functions that I commonly use
#libraries

library(viper)
library(gplots)
library(DESeq2)
library(randomForest)
library(varSelRF)
library(rfUtilities)
library(survival)
library(rms)
library(survminer)
library(SummarizedExperiment)
library(Biobase)
library(ggplot2)
library(cluster) #for running pam()
library(limma)
library(fmsb)
library(ROCR)

library(MASS)
library(glmnet)
library(caret)
library(pROC)

library(tsne)
library(atools)
library(clusterpam)
library(MKmisc)
library(ConsensusClusterPlus)

library(DPpackage)

library(metaseqR)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(rgl) #for plot3d
#####Package installation##
#If running a new installation of R
install.packages('readxl')
install.packages('rgl')
install.packages('DPpackage')

source("http://bioconductor.org/biocLite.R")

biocLite('viper')
biocLite('gplots')
biocLite('DESeq2')
biocLite('randomForest')
biocLite('varSelRF')
biocLite('rfUtilities')
biocLite('survival')
biocLite('rms')
biocLite('survminer')
biocLite('SummarizedExperiment')
biocLite('Biobase')

biocLite("ggplot2")
biocLite("cluster")
biocLite("limma")
biocLite("fmsb")
biocLite("ROCR")


biocLite("MASS")
biocLite("glmnet")
biocLite("caret")
biocLite("pROC")

biocLite("tsne")
biocLite("atools")
biocLite("clusterpam")
biocLite("MKmisc")
biocLite('ConsensusClusterPlus')

biocLite('metaseqR')

biocLite("org.Hs.eg.db")
biocLite('GenomicFeatures')
biocLite('TxDb.Hsapiens.UCSC.hg38.knownGene')
#############################


# source('./hallmarks.radar.plot.R')



#functions:
# 1- PCA (runs PCA on either all the genes of a matrix or a subset of it, then shows the output as a graph)
# 2- Heatmap (uses heatmap.2 package to draw a heatmap with an optional colors_hm parameter to set the column bar columns)
# 3- Dictionary (to convert gene names between 'genenames', 'entrez' and 'ensemble' formats)
PCA <- function (matrix,title='PCA', pt_colors='black', filter=0, show_scree=FALSE, top_genes = nrow(matrix), component_1 = 1, component_2 = 2, lines = FALSE, output = NULL, pch = 16, show_legend = T, return_PCA_matrix = F, return_full_PCA_object = T){
  if (typeof(matrix) != 'matrix'){
    matrix <- as.matrix(matrix)
    total_nb_genes <- nrow(matrix)
  }
  if (filter == 0 & top_genes == nrow(matrix)){
    print(paste0('PCA on 100% of all genes (',total_nb_genes,' genes).'))
    PCA <- prcomp(t(matrix))
  }
  if (filter !=0 & top_genes == nrow(matrix)){
    selected_nb_genes <- nrow(matrix[rowVars(matrix) > filter*max(rowVars(matrix)) , ])
    print(paste0('PCA on ', round(100*selected_nb_genes/total_nb_genes, 1) ,'% of all genes (',selected_nb_genes,' of ', total_nb_genes,' genes).'))
    PCA <- prcomp(t(matrix[rowVars(matrix)>filter*max(rowVars(matrix)),]))
  }
  if (top_genes < nrow(matrix)){
    selected_nb_genes <- top_genes
    print(paste0('PCA on ', round(100*selected_nb_genes/total_nb_genes, 1) ,'% of all genes (',selected_nb_genes,' of ', total_nb_genes,' genes).'))
    matrix <- matrix[order(rowVars(matrix),decreasing = TRUE)[1:top_genes],]
    PCA <- prcomp(t(matrix))
  }
  if (top_genes > nrow(matrix)){
    #print('Error: you selected more genes than there are in the matrix!')
    stop('You selected more genes than there are in the matrix!\nRerun the command with a smaller number of genes under as top_genes condition, or leave blank. ')
  }
  
  
  if (show_scree == F){
    #par(mfrow=c(1, 1))
    plot(x = PCA$x[ , component_1], y=PCA$x[ , component_2], col=pt_colors, pch=pch, xlab=paste0('PC',component_1,' - Variance explained: ',round(100*(PCA$sdev[component_1]/sum(PCA$sdev)),2),'%'), ylab=paste0('PC',component_2,' - Variance explained: ',round(100*(PCA$sdev[component_2]/sum(PCA$sdev)),2),'%'), main = title)
    
    if (show_legend){
      legend(x = 'bottomright', legend = levels(as.factor(pt_colors)), fill = levels(as.factor(pt_colors)))  
    }
    
    if (lines == T){
      for (i in 1:35){
        segments(x0 = PCA$x[i,1], x1 = PCA$x[(i+35),1], y0 = PCA$x[i,2], y1 = PCA$x[(i+35),2])#, col = pt_colors)
      }
    }
  }
  
  
  
  if (show_scree == T){
    par(mfrow=c(1, 2))
    plot(x = PCA$x[ , component_1], y=PCA$x[ , component_2], col=pt_colors, pch=pch, xlab=paste0('PC',component_1,' - Variance explained: ',round(100*(PCA$sdev[component_1]/sum(PCA$sdev)),2),'%'), ylab=paste0('PC', component_2 ,' - Variance explained: ',round(100*(PCA$sdev[component_2]/sum(PCA$sdev)),2),'%'), main = title)
    PCA$sdev <- sapply(PCA$sdev, function(x) (100*x)/sum(PCA$sdev))
    
    if (show_legend){
      legend(x = 'bottomright', legend = levels(as.factor(pt_colors)), fill = levels(as.factor(pt_colors)))  
    }
    
    barplot(height = PCA$sdev[1:10], main = 'Scree plot for PCA', ylab = '% Variance explained', names.arg = sprintf("PC%d", 1:10))
  }
  
  if (is.null(output) == FALSE){
    pdf(file = output, width = 9.375, height = 9.375, title = title)
    dev.off()
  }
  
  if (return_PCA_matrix == T){
    PCA_res = list()
    PCA_res[['x']] = PCA$x[ , component_1]
    PCA_res[['y']] = PCA$x[ , component_2]
    return(PCA_res)
  }
  
  if (return_full_PCA_object == T){
    return(PCA)
    
  }
  
}

MDS <- function(input,  plot_title = 'MDS', pt_colors = 'black', metric = T){
  input = t(input)
  
  if (metric == T){
    d = dist(input)
    fit = cmdscale(d, eig = TRUE, k = 2) 
    
    plot(x = fit$points[,1], y = fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",
         main = plot_title, col = pt_colors)  
  }
  
  if (metric == F){
    d = dist(input) # euclidean distances between the rows
    fit = isoMDS(d, k=2) # k is the number of dim
    
    plot(x = fit$points[,1], y = y <- fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",
         main = plot_title, col = pt_colors)
  }
  
}

HM <- function(input, subset = 1:ncol(input), nn= NULL,
               colors_hm =NULL, colors_hm_y = NULL, scale_color = NULL, scale_row_colors = F, 
               vipersimilarity = FALSE, filter = 1, filter_top = TRUE, 
               input_hc = NULL, return_hc = TRUE, return_hr = FALSE, return_hm = FALSE,
               cexRow = 1, cexCol = 0.3, key = FALSE)
  {
  #filter: value from 0 to 1, 1 will keep all genes, 0.01 will keep the 1% top genes sorted by variance
  #vipersimilarity: TRUE will compute the vipersimilarity matrix and cluster accordingly
  
    
  if (filter_top == TRUE){ #will select the top genes (the genes with the most variance)
    input <- input[order(rowVars(input),decreasing = T)[1:(nrow(input)*filter)],]  
  }
  
  if (filter_top == FALSE){ #will select the bottom genes (the genes with the least variance)
    input <- input[order(rowVars(input),decreasing = F)[1:(nrow(input)*filter)],]
  }
  
  
  
  hc <- hclust(as.dist(1-cor((input[,subset]), method='spearman')), method ='complete')
  hr <- hclust(as.dist(1-cor(t(input), method='spearman')), method ='complete')
  
  if (vipersimilarity != FALSE){
    if (vipersimilarity == TRUE){
      vipersimilarity = 'greater' #defaults to 'greater' if user picks 'TRUE'. other options: 'less' and 'two.sided'.
    }
    vps <- viperSimilarity(x = input, method = vipersimilarity, nn = nn)
    class(vps) <- 'matrix'
    vps <- as.data.frame(vps)
    if (is.null(input_hc == T)) {hc = hclust(as.dist(1-cor((vps[,subset]), method='spearman')),method='complete')}
    if (is.null(input_hc == F)) {hc = input_hc}
  }
  
  
  
  
  if (scale_row_colors == T){
    input = scale(x = t(input), center = T, scale = T)
    input = t(input)
  }
  if (return_hm == TRUE){
    return_hc = F  #set only one condition to be true at a time. return_hm will return the heatmap as an object
  }
  
  
  
  if (return_hm == FALSE & return_hc == FALSE & is.null(colors_hm) == F){
    hm = heatmap.2(x = input[,subset], col = colorRampPalette(c('blue','white','red'))(n=20), breaks = scale_color, tracecol = NULL, dendrogram = 'none', Rowv = F, Colv = F, ColSideColors = as.character(as.factor(colors_hm)), cexRow = cexRow, cexCol = cexCol, key = key)
    return(hm)
  }
  
  if (return_hm == FALSE & return_hc == FALSE & is.null(colors_hm) == T){
    hm = heatmap.2(x = input[,subset], col = colorRampPalette(c('blue','white','red'))(n=20), tracecol = NULL, dendrogram = 'none', Rowv = F, Colv = F, cexRow = cexRow, cexCol = cexCol, key = key)
    return(hm)
  }
  
  
  
  #scale color in this form: scale_color = seq(-3,3,length.out = 21)
  if (is.null(colors_hm) == TRUE){
    hm = heatmap.2(input[,subset], col = colorRampPalette(c('blue','white','red'))(n=20), breaks = scale_color, tracecol = NULL, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), cexRow = cexRow, cexCol = cexCol, key = key)
  }
  
  print(ncol(input[,subset]))
  
  if (is.null(colors_hm) == FALSE & is.null(colors_hm_y) == FALSE){
    hm = heatmap.2(input[,subset], col = colorRampPalette(c('blue','white','red'))(n=20), breaks = scale_color, tracecol = NULL, ColSideColors = as.character(as.factor(colors_hm)), RowSideColors = as.character(as.factor(colors_hm_y)) ,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), cexRow = cexRow, cexCol = cexCol, key = key)  
  }
  
  if (is.null(colors_hm) == FALSE & is.null(colors_hm_y) == TRUE){
    hm = heatmap.2(input[,subset], col = colorRampPalette(c('blue','white','red'))(n=20), breaks = scale_color, tracecol = NULL, ColSideColors = as.character(as.factor(colors_hm)), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), cexRow = cexRow, cexCol = cexCol, key = key)  
  }

  if (return_hr == TRUE){
    return(hr)
  }  
  if (return_hm == TRUE){
    return(hm)
  }
  if (return_hc == TRUE){
    return(hc)
  }

}

rank_rank_transformation <- function(GEM){
  col_ranks = apply(GEM, 2, rank)
  row_medians = apply(col_ranks, 1, median)
  row_mads = apply(row_medians, 1, mad)
  return( (col_ranks - row_medians)/row_mads )
}

double_rank_transformation <- function(tpm){
  rank <- apply(tpm, 2, rank) ; print(paste0('rank1: ', rank[32516,]))
  median <- apply(rank, 1, median) ; print(paste0('median: ', median[32516] ))
  mad <- apply(rank, 1, mad) ; print(paste0('mad: ', mad[32516]))
  rank <- (rank - median)/mad ; #print(paste0('rank2: ', rank[32516] ))
  return(rank)
}



Genename_Conv <- function(source = 'entrez', target = 'genenames', vector_to_translate){
  if (exists('dict') == F){
    print('dictionary not in ls()')
    assign(x = 'dict', value = '/Users/alexandre/Desktop/Columbia/CalifanoLab/scripts/R/Genename_dictionary.txt', envir = .GlobalEnv)
    assign(x = 'dict', value = read.delim(dict, sep='\t', header=FALSE, stringsAsFactors = FALSE), envir = .GlobalEnv)
    genenames <<- dict[,1]
    entrez <<- dict[,2]
    ensembl <<- dict[,3]
    ensembl[which(ensembl == '')] <- NA
  }
  
  if ((source != 'entrez' & source != 'genenames' & source != 'ensembl') | (target != 'entrez' & target != 'genenames' & target != 'ensembl')) {
    stop('Error: source and target need to be either "entrez", "genenames" or "ensembl".')
  }
  
  if (vector_to_translate[1] %in% get(source) == FALSE){
    stop('Error: the vector to translate is not in the right format (you may be trying to translate an entrez vector by the gene names are in another format)')
  }
  
  return( get(target)[match(vector_to_translate, get(source))]  )
  
}

convert_ensembl_to_entrez_matrix <- function(matrix_with_ensembl_id){
  vector_row_names = row.names(matrix_with_ensembl_id)
  vector_row_names = Genename_Conv('ensembl', 'entrez', vector_row_names)
  list_NAs = which(is.na(vector_row_names))
  list_empties = which(vector_row_names == '')
  list_to_remove = c(list_NAs, list_empties)
  mat_res = matrix_with_ensembl_id
  
  mat_res = mat_res[setdiff(1:nrow(mat_res), list_to_remove),] ; print(dim(mat_res))
  vector_row_names = vector_row_names[setdiff(1:length(vector_row_names), list_to_remove)]
  row.names(mat_res) = vector_row_names
  return(mat_res)
}

convert_genenames_matrix <- function(original_matrix, source_id = 'ensembl', dest_id = 'entrez'){
  vector_row_names = row.names(original_matrix)
  vector_row_names = Genename_Conv(source_id, dest_id, vector_row_names)
  list_NAs = which(is.na(vector_row_names))
  list_empties = which(vector_row_names == '')
  list_to_remove = c(list_NAs, list_empties)
  mat_res = original_matrix
  
  mat_res = mat_res[setdiff(1:nrow(mat_res), list_to_remove),] ; print(dim(mat_res))
  vector_row_names = vector_row_names[setdiff(1:length(vector_row_names), list_to_remove)]
  row.names(mat_res) = vector_row_names
  return(mat_res)
}





Translate_Marina <- function(msviper_object, source_ID = 'entrez', target_ID = 'genenames'){
  if (class(msviper_object) != 'msviper'){
    print('Error: the msviper_object has to be an object of class "msviper"')
    break
  }
  names(msviper_object$es$nes) <- Genename_Conv(source_ID, target_ID, names(msviper_object$es$nes))
  names(msviper_object$es$p.value) <- Genename_Conv(source_ID, target_ID, names(msviper_object$es$p.value))
  names(msviper_object$es$size) <- Genename_Conv(source_ID, target_ID, names(msviper_object$es$size))
  names(msviper_object$regulon) <- Genename_Conv(source_ID, target_ID, names(msviper_object$regulon))
  row.names(msviper_object$signature) <- Genename_Conv(source_ID, target_ID, row.names(msviper_object$signature))
  row.names(msviper_object$es$nes.bt) <- Genename_Conv(source_ID, target_ID, row.names(msviper_object$es$nes.bt))
  if (source_ID == 'entrez'){
    print('Translating gene names in msVIPER regulon object from entrez ID to genenames...')  
  }
  if (source_ID == 'genenames'){
    print('Translating gene names in msVIPER regulon object from Genenames to entrez ID')
  }
  
  pb <- txtProgressBar(max=length(names(msviper_object$es$nes)), style=3)
  x <- 0
  for (i in names(msviper_object$es$nes)){
    x = x + 1
    setTxtProgressBar(pb,x)
    names(msviper_object['regulon'][[1]][[i]]$tfmode) <- Genename_Conv(source = source_ID, target = target_ID, vector_to_translate = names(msviper_object['regulon'][[1]][[i]]$tfmode))
  }
  
  return(msviper_object)
}

k_fold_randomforest <- function(k, input_matrix, sample_labels, ntree = 1000){
  
  nb_samples <- length(sample_labels)
  nb_samples_in_k <- floor(nb_samples / k)
  
  if (nb_samples_in_k < 1){
    print('Error: k too large')
    return()
  }
  
  pool_available_samples <- sample_labels
  for (i in 1:k){
    
    assign(x = paste0('test_set_',i), value = sample(x = pool_available_samples, size = nb_samples_in_k, replace = F))
    pool_available_samples <- pool_available_samples[setdiff(names(pool_available_samples), names(get(paste0('test_set_',i))))]
    
    if (i == k & (k*nb_samples_in_k) < nb_samples){
      assign(x = paste0('test_set_',i), value = c(get(paste0('test_set_',i)), pool_available_samples))
    }
  assign(x = paste0('training_set_',i), value = sample_labels[setdiff(names(sample_labels) , names(get(paste0('test_set_',i))))])

  }
 
  
  for (i in 1:k){
    assign(x = paste0('rf_',i),
           value = randomForest(x = input_matrix[names(get(paste0('training_set_',i))),],
                                y = as.factor(get(paste0('training_set_',i))),
                                xtest = input_matrix[names(get(paste0('test_set_',i))),],
                                ytest = as.factor(get(paste0('test_set_',i))),
                                ntree = ntree, do.trace = T))
  }
  
  
   
}

clinical_data_check <- function(Clinical_Data_Matrix, Name_Column, Patient_Names){
  for (i in 1:ncol(Clinical_Data_Matrix)){
    if (typeof(Clinical_Data_Matrix[,i]) == 'integer'){
      print(colnames(Clinical_Data_Matrix)[i])
      print(summary(Clinical[,i][Clinical_Data_Matrix[,Name_Column] %in% Patient_Names]))
      print('-----------------------------------------------------------------------------------------------')
    }
    if (typeof(Clinical_Data_Matrix[,i]) == 'character'){
      print(colnames(Clinical_Data_Matrix)[i])
      print(table(Clinical[,i][Clinical_Data_Matrix[,Name_Column] %in% Patient_Names]))
      print('-----------------------------------------------------------------------------------------------')
    }
  }
}
#Example:
# clinical_data_check(Clinical_Data_Matrix = Clinical, Name_Column = 1, Patient_Names = colnames(vp_Diag)[which(vp_Diag['DOT1L',] < (-5))])


#ROC curve to use with random forest objects
ROC_curve <- function(rf, labels, plot_title = 'ROC plot'){
  library(ROCR)
  predictions=as.vector(rf$votes[,2])
  pred=prediction(predictions,labels)
  perf_AUC=performance(pred,"auc") 
  AUC=perf_AUC@y.values[[1]]
  perf_ROC=performance(pred,"tpr","fpr")
#  plot(perf_ROC, main= plot_title)
#  abline(a = 0, b = 1)
#  text(x = 0.8, y = 0.1, labels = paste0('AUC = ', round(AUC,2)))
  return(AUC)
}

RF <- function(matrix, subset1, subset2, subset_name_1 = 'subset1', subset_name_2 = 'subset2', ntree = 500, variable_selection = F, ROC_curve_title = 'ROC plot'){
  matrix <- t(matrix[,c(subset1,subset2)])
  nb_patients_subset_1 <- length(subset1)
  nb_patients_subset_2 <- length(subset2)
  classification_labels <- as.factor(
    c(
      rep(x = subset_name_1, times = nb_patients_subset_1),
      rep(x = subset_name_2, times = nb_patients_subset_2)
    )
  )
  
  if (variable_selection == F){
    output <- randomForest(x = matrix, y = classification_labels, ntree = ntree, do.trace = T)  
  }
  
  if (variable_selection == T){
    output <- varSelRF(xdata = matrix, Class = classification_labels, ntree = ntree ,verbose = T)
    return(output)
  }
  
  ROC_curve(rf = output, labels = classification_labels, plot_title = ROC_curve_title)
  return(output)
  
}


######################################################################################
#Function to pick the best individual features by AIC using logistic regression with all the features individually
# matrix: input matrix (viper or gene expression matrix, patients on the columns and genes/proteins on the rows)
# subtype: must enter 'M2', 'M4' or 'M5'
# nb_features: the length of the output, will show the top hits (default is 10)
# show_heatmap: will show a heatmap of the top hits 
#will run logistic regression with all features one by one and return a list of AICs that can be sorted. The top AICs are the features that classify best
glm_subtype_logit_class <- function(matrix, subtype, nb_features = 10, show_heatmap = T, show_top_features = F){
  NR_length = length(non.relapse[FAB_NR == subtype])
  R_length = length(relapse[FAB_R == subtype])
  vp_subtype = rbind(matrix, 'relapse_status' = c(  rep(0, NR_length), rep(1, R_length)  ))
  glm_vp_subtype = sapply(X = 1:nrow(matrix), FUN = function(i) (glm(vp_subtype[nrow(vp_subtype),] ~ vp_subtype[i,], family = binomial(link = 'logit')))$deviance  )
  
  features = row.names(matrix)[order(glm_vp_subtype, decreasing = F)[1:nb_features]]
  
  if(show_heatmap == T){
    matrix = matrix[features,]
    matrix = t(matrix)
    matrix = scale(x = matrix, center = T, scale = T)
    matrix = t(matrix)
    HM(input = matrix, colors_hm = c(  rep('blue', NR_length), rep('red', R_length)  ), return_hc = F)  
  }
  
  if (show_top_features == T){
    par(mfrow=c(3,3))
    for (i in 1:9){
      plot(x = vp_subtype[features[i],], y = vp_subtype[nrow(vp_subtype),], main = features[i], xlab = 'VIPER NES', ylab = 'Relapse(1) vs non-relapse (0)')
      abline(v = 0, col = 'red')
    }  
  }
  
  print(features)
  par(mfrow=c(1,1))
  return (features)  
}


######################################################################################
######################################################################################
# Cross-validation LASSO logistic regresssion function to be run on VIPER matrix to separate relapse and non-relapse patients

# input: VIPER expression matrix, patients in columns, proteins in rows
# response_list: 1 x n vector of responses. Must be binary (0 or 1).

# alpha: value for alpha (1 for lasso, 0 for ridge regression, 0.5 for elastic net) Can pick any number between 0 and 1.
# optimize_alpha: will run cv.glmnet with different values of alpha and will plot mean-squared-error vs log(lambda)
# alpha_range: a list of alpha values to supply for optimizing. Can be in this format alpha_range = c(0.1,0.2,0.3,0.4) or a seq(from = 0.1, to = 1, length.out = 9)

# LOOCV: boolean, leave-one-out-cross-validation. If TRUE will run LOOCV instead of k-fold CV.
# k: number of folds for cross-validation, if LOOCV is FALSE. For example k = 5 will run 5-fold cross-validation

# repeat_cvglmnet: integer (n). Will repeat cv.glmnet n times using the conditions defined below (k-fold cv or LOOCV, alpha etc) Will return a list of features with the number of times they show up in the repeated runs
# draw_ROC_repeat: Will draw an ROC curve. Each run will have its own curve.
# draw_ROC_repeat_null: Will draw 1000 ROC curves with the labels shuffled (like null_model_ROC below)

# show_HM: will show a heatmap of the top genes/features found by logistic regression and all the patients. A color bar on top will indicate the response
# centered_NESes: will scale each row of the heatmap so that one can more easily visualize the difference in activity between the two responses

# graph_title: this title will appear on the ROC graph
# show_ROC: will show the receiver-operator curve for the prediction
# null_model_ROC: will display null-model ROCs with 1000 permutations of the response labels

# return_model: will return the glmnet model. You can save it in a variable and then obtain the beta coefficients and other features of the model
# monte_carlo: will run a monte-carlo cross-validation instead of k-fold cross-validation (if "NULL" will run k-fold CV, if not NULL, will run MCCV with the desired proportion used as test set (input from 0 to 1))


cv.glmnet.FAB.AML <- function(input, response_list, single_predictor = F,
                              alpha = 1, optimize_alpha = F, alpha_range = ..., 
                              LOOCV = T, k = 5, monte_carlo = FALSE,
                              repeat_cvglmnet = 1, draw_ROC_repeat = F, draw_ROC_repeat_null = F,
                              show_HM = F, centered_NESes = T, 
                              graph_title = 'ROC', show_ROC = T, null_model_ROC = T, 
                              return_model = F){
  input = t(input)
  input = as.matrix(input)
  response_list <- as.factor(response_list)
  lambda_list = list()
  
  if (single_predictor == TRUE){
    glm_single_predictor = sapply(X = 1:nrow(input), FUN = function(i) (glm(response_list ~ input[i,], family = binomial(link = 'logit')))$aic  )
    return(glm_single_predictor)
    #coefficients.lasso = glm_single_predictor
  }
  
  if (optimize_alpha == T){
    par(mfrow = rep(ceiling(sqrt(length(alpha_range))), 2) )
    for (i in 1:length(alpha_range)){
      
      if (LOOCV == F){
        cvfit = cv.glmnet(x = as.matrix(input), y = response_list, type.measure = "class", nfolds = k, alpha = alpha_range[i], grouped = F, family = 'binomial')  
      }
      if (LOOCV == T){
        cvfit = cv.glmnet(x = as.matrix(input), y = response_list, type.measure = "class", nfolds = length(response_list), alpha = alpha_range[i], grouped = F, family = 'binomial')    
      }
      plot(cvfit, main = paste0('alpha = ', alpha_range[i])  )
      abline(h = cvfit$cvm[which(cvfit$lambda == cvfit$lambda.min)] , col = 'green')
      abline(h = 0.5 , col = 'red')
      lambda_list[[as.character(alpha_range[i])]] <- cvfit$lambda.min
    }
    par(mfrow = c(1,1))
    return(lambda_list)
  }
  
  if (repeat_cvglmnet > 1){
    print(paste0('Repeating cv.glmnet ', as.character(repeat_cvglmnet), ' times with ', as.character(k), '-fold cross-validation.'))
    pb = txtProgressBar(min = 0, max = repeat_cvglmnet, initial = 1, style = 3)
    list_coeff_runs <- list()
    for (i in 1:repeat_cvglmnet){
      setTxtProgressBar(pb = pb, value = i)
      
      cvfit = cv.glmnet(x = as.matrix(input), y = response_list, type.measure = "class", nfolds = k, alpha = alpha, grouped = F, family = 'binomial')
      lasso.fit = cvfit$glmnet.fit
      lambda.min = cvfit$lambda.min
      beta.matrix.index = which(lasso.fit$lambda == lambda.min)
      coefficients.lasso = lasso.fit$beta[,beta.matrix.index][lasso.fit$beta[,beta.matrix.index] != 0 ]
      list_coeff_runs[[i]] = coefficients.lasso
      
      if (draw_ROC_repeat == TRUE){
        if (i == 1){
          pred_test = predict(object = cvfit, newx = as.matrix(input), s = 'lambda.min', type = 'response')
          pred <- prediction(pred_test, response_list)
          perf <- performance(pred,"tpr","fpr")
          performance(pred,"auc")
          plot(perf, col="black", main = graph_title)
          abline(a = 0, b = 1, col = 'red')  
        }
        if (i > 1){
          pred_test = predict(object = cvfit, newx = as.matrix(input), s = 'lambda.min', type = 'response')
          pred <- prediction(pred_test, response_list)
          perf <- performance(pred,"tpr","fpr")
          lines(x = perf@x.values[[1]], y = perf@y.values[[1]], col = 'black'  )
        }
      }
      
    }
    
    if (draw_ROC_repeat_null == TRUE){
      random_labels_shuffling = matrix(nrow = 1000, ncol = length(response_list))
      for (i in 1:1000){
        random_labels_shuffling[i,] = sample(x = response_list, size = length(response_list), replace = T)
      }  
      
      for (i in 1:1000){
        pred_r <- prediction(pred_test, random_labels_shuffling[i,])
        perf_r <- performance(pred_r,"tpr","fpr")
        lines(x = perf_r@x.values[[1]], y = perf_r@y.values[[1]], col = 'grey', lty = 3  )
      }
      # lines(x = perf@x.values[[1]], y = perf@y.values[[1]], col="black")
      abline(a = 0, b = 1, col = 'red')
    }
    
    
    res = table(names(unlist(list_coeff_runs)))
    return(res[order(res, decreasing = T)])
  }
  
  #If you supply a value for alpha and set optimize_alpha = FALSE, code starts here:
  if (LOOCV == F){
    cvfit = cv.glmnet(x = as.matrix(input), y = response_list, type.measure = "class", nfolds = k, alpha = alpha, grouped = F, family = 'binomial')
    cvfit$cvm[which(cvfit$lambda == cvfit$lambda.min)]
  }
  if (LOOCV == T){
    cvfit = cv.glmnet(x = as.matrix(input), y = response_list, type.measure = "class", nfolds = length(response_list), alpha = alpha, grouped = F, family = 'binomial')    
  }
  
  ######################################################################################################################################################################################################
  ######################################################################################################################################################################################################
  
  if (monte_carlo == TRUE){
    #method : LGOCV = leave group out CV = Monte Carlo CV
    #number: number of iterations
    #p: the training percentage for LGOCV
    trControl <- trainControl(method = "LGOCV", number = 1000, p = 0.75, verboseIter = T, savePredictions = 'final', classProbs = TRUE)
    glmnet.obj <- train(y = as.factor(make.names(response_list)), x = input, method = 'glmnet', family = 'binomial', trControl = trControl)    
    
    return(row.names(glmnet.obj$modelInfo$varImp(object = glmnet.obj$finalModel, lambda = glmnet.obj$finalModel$lambdaOpt))[which(glmnet.obj$modelInfo$varImp(object = glmnet.obj$finalModel, lambda = glmnet.obj$finalModel$lambdaOpt)!=0)])
    
  }
  
  #############################################################################################################################################################################################
  #############################################################################################################################################################################################
  
  plot(cvfit)
  lasso.fit = cvfit$glmnet.fit
  lambda.min = cvfit$lambda.min
  beta.matrix.index = which(lasso.fit$lambda == lambda.min)
  coefficients.lasso = lasso.fit$beta[,beta.matrix.index][lasso.fit$beta[,beta.matrix.index] != 0 ]
  print(coefficients.lasso)
  coefficients.lasso = names(coefficients.lasso)
  print(coefficients.lasso)
  
  
  if (show_HM == TRUE){
    
    if (centered_NESes == T){
      input_HM = scale(x = input[,coefficients.lasso], center = T, scale = T)
      input_HM = t(input_HM)
    }
    
    if (centered_NESes == F){
      input_HM = t(input)  
    }    
    
    colors_HM = NULL
    colors_HM[which(response_list == 0)] <- 'blue'
    colors_HM[which(response_list == 1)] <- 'red'
    
    HM(input = input_HM[coefficients.lasso,c(which(colors_HM == 'blue'), which(colors_HM == 'red'))], colors_hm = colors_HM[c(which(colors_HM == 'blue'), which(colors_HM == 'red'))], return_hc = F)  
  }
  
  #ROC curve + null model
  if (show_ROC == TRUE){
    pred_test = predict(object = cvfit, newx = as.matrix(input), s = 'lambda.min', type = 'response')
    pred <- prediction(pred_test, response_list)
    perf <- performance(pred,"tpr","fpr")
    performance(pred,"auc")
    plot(perf, col="black", main = graph_title)
    abline(a = 0, b = 1, col = 'red')  
  }
  
  if (null_model_ROC == TRUE){
    random_labels_shuffling = matrix(nrow = 1000, ncol = length(response_list))
    for (i in 1:1000){
      random_labels_shuffling[i,] = sample(x = response_list, size = length(response_list), replace = T)
    }  
    
    for (i in 1:1000){
      pred_r <- prediction(pred_test, random_labels_shuffling[i,])
      perf_r <- performance(pred_r,"tpr","fpr")
      lines(x = perf_r@x.values[[1]], y = perf_r@y.values[[1]], col = 'grey', lty = 3  )
    }
    lines(x = perf@x.values[[1]], y = perf@y.values[[1]], col="black")
    abline(a = 0, b = 1, col = 'red')
  }
  
  if (return_model == T){
    return(cvfit)
  }
  
  
  return(coefficients.lasso)  
  
  
}





ROC_glmnet_monte_carlo <- function(caret_object, viper_matrix, response_list, graph_title = 'ROC curve', reverse_pred = F, null_model = 1000, show_ROC = T, return_AUC = F){
  
  if (reverse_pred == TRUE){
    predictions = caret_object$pred$X0  
  }
  
  if (reverse_pred == FALSE){
    predictions = caret_object$pred$X1
  }
  
  classification_labels <- caret_object$pred$rowIndex
  classification_labels[classification_labels %in% which(response_list == names(table(response_list)[1]))] <- 'condition_1'
  classification_labels[classification_labels %in% which(response_list == names(table(response_list)[2]))] <- 'condition_2'
  
  classification_labels = make.names(classification_labels)
  #print(classification_labels)
  pred <- prediction(predictions = predictions, classification_labels)
  perf_AUC=performance(pred,"auc") 
  AUC=perf_AUC@y.values[[1]]
  perf_ROC=performance(pred,"tpr","fpr")
  
  if (show_ROC == T){
    plot(perf_ROC, main = graph_title)
    text(x = 0.8, y = 0.2, paste0('AUC = ', round(AUC, 2)))
    
    distribution_random_AUCs <- list()
    pb <- txtProgressBar(min = 1, max = null_model, style = 3)
    for (i in 1:null_model){
      setTxtProgressBar(pb, i)
      random_classification <- sample(x = c('bad', 'good'), size = length(classification_labels), replace = T)
      pred_random <- prediction(predictions, random_classification)
      perf_random = performance(pred_random, 'tpr', 'fpr')
      random_points=sample(x = 1:length(predictions), size = 1, replace = F)
      decision_to_plot_point = sample(x = 1:100, size = 1)
      if (decision_to_plot_point < 101){
        points(x = perf_random@x.values[[1]][random_points], y = perf_random@y.values[[1]][random_points], col = 'grey')    
      }
      
      distribution_random_AUCs[i] <- performance(pred_random, 'auc')@y.values[[1]]
    }
    abline(a = 0, b = 1, col = 'red')
    pvalue = 1-pnorm((AUC - mean(unlist(distribution_random_AUCs)))/sd(unlist(distribution_random_AUCs)))
    if (pvalue < 1/1000){
      pvalue = 1/1000
    }
    text(x = 0.8, y = 0.1, paste0('p-value = ', pvalue ) )
  }

  
  if (return_AUC == T){
    return(AUC)  
  }
  
}

monte_carlo_cross_validation_lasso <- function(input, response_list, alpha = 1, lambda_1SE = TRUE, show_ROC = F, return_best_features = F, return_full_model = T, lambda_value = NULL, compute_best_lambda = 'kfold', verbose = T){
  trControl <- trainControl(method = "LGOCV", number = 1000, p = 0.75, verboseIter = T, savePredictions = TRUE, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  input = t(input)
  list_lambdas <- list()
  
  if (is.null(lambda_value) == T){
    
    if (verbose == T){
      print('Computing optimal lambda for MCCV (lambda.1se)')  
    }
    
    if (compute_best_lambda == 'kfold'){
      pb = txtProgressBar(min = 1, max = 20, style = 3)
      for (i in 1:20){
        if (verbose == T)
        {
          setTxtProgressBar(pb, i)
        }
        
        lambda_max <- cv.glmnet(x = input, y = as.factor(make.names(response_list)), family = 'binomial', alpha = alpha, nfolds = 5, type.measure = 'class')
        if (lambda_1SE == T){
          list_lambdas[[i]] <- lambda_max$lambda.1se
        }
        if (lambda_1SE == F){
          list_lambdas[[i]] <- lambda_max$lambda.min
        }
      }
    }
    
    
    #Use LOOCV if you have a small sample size. Make sure you have at least 2 of each class, otherwise some folds will end up with 0 of the class of interest and will return an error
    if (compute_best_lambda == 'LOOCV'){
      lambda_max <- cv.glmnet(x = input, y = as.factor(make.names(response_list)), family = 'binomial', alpha = alpha, nfolds = length(response_list), type.measure = 'class')
      if (lambda_1SE == T){
        list_lambdas[[1]] <- lambda_max$lambda.1se  
      }
      if (lambda_1SE == F){
        list_lambdas[[1]] <- lambda_max$lambda.min  
      }
      
    }
    
    
    if (verbose == T){
      if (lambda_1SE == T){
        print(paste0('lambda.1se = ', median(unlist(list_lambdas)))  )
      }
      
      if (lambda_1SE == F){
        print(paste0('lambda.min = ', median(unlist(list_lambdas)))  )
      }  
    }
    
    
  }
  
  if (is.null(lambda_value) == F){
    list_lambdas = rep(lambda_value, 10)
  }
  
  print('Computing the MCCV model')
  
  MCCV_model <- train(y = make.names(as.factor(response_list)), 
                      x = input, 
                      method = 'glmnet', family = 'binomial', metric = 'ROC',
                      trControl = trControl,
                      tuneGrid = expand.grid(.alpha = alpha, .lambda = median(unlist(list_lambdas))  ))
  
  if (show_ROC){
    ROC_glmnet_monte_carlo(caret_object = MCCV_model, viper_matrix = t(input), response_list = response_list, null_model = 10)
  }
  
  if (return_best_features == T & return_full_model == T){
    MCCV_model$top_features = row.names(MCCV_model$modelInfo$varImp(object = MCCV_model$finalModel, lambda = MCCV_model$finalModel$lambdaOpt))[which(MCCV_model$modelInfo$varImp(object = MCCV_model$finalModel, lambda = MCCV_model$finalModel$lambdaOpt)!=0)]
    return(MCCV_model)
  }
  
  if (return_best_features){
    return(row.names(MCCV_model$modelInfo$varImp(object = MCCV_model$finalModel, lambda = MCCV_model$finalModel$lambdaOpt))[which(MCCV_model$modelInfo$varImp(object = MCCV_model$finalModel, lambda = MCCV_model$finalModel$lambdaOpt)!=0)])
  }
  
  if (return_full_model){
    return(MCCV_model)
  }
}

CC <- function(input, distance_metric = 'spearman', cluster_alg = 'pam', prop_features = 1, results_path = 'ConsensusClustering/', show_silhouette = T){
  consensus_res = ConsensusClusterPlus(d = as.matrix(input), maxK = 10, reps = 1000, distance = distance_metric, clusterAlg = cluster_alg, verbose = T, plot = 'png', writeTable = F, title = results_path, pFeature = prop_features)
  if (show_silhouette){
    par(mfrow = c(3,3))
    for (i in 2:10){
      plot(silhouette(x = consensus_res[[i]]$consensusClass, dist = as.dist(1-cor(as.matrix(input), method = distance_metric))), main = '')
      print ( mean(silhouette(x = as.numeric(consensus_res[[i]]$consensusClass), dist = as.dist(1 -cor(x = as.matrix(input), method = distance_metric)))  [, 'sil_width'] ) )
      #print(mean(silhouette(x = consensus_res[[i]]$consensusClass, dist = as.dist(1-cor(as.matrix(input), method = distance_metric)))))
    }
  }
  return(consensus_res)
}

Plot_Clinical_Features <- function(clinical_features_matrix = NULL, clin_features = 1:nrow(clinical_features_matrix), cluster_assignments, arbitrary_clusters = NULL, p_val_threshold = 0.05, FDR = T){
  
  nb_features = length(clin_features) #total number of features to look at
  nb_clusters = length(table(cluster_assignments)) #total number of clusters in the cluster_assignment object
  
  #Formation of clusters and of the clinical data matrix
  #Assigns each patient to cluster 1, 2, 3, ...
  #Assemble the clinical features (the rows selected in clin_features) from the matrix clinical_features_matrix
  #End result: a matrix with patient names ordered by cluster (columns) and clinical features (rows)
  COUNTER = 0
  for (cluster_id in (1:nb_clusters)){
    if (COUNTER == 0){
      patient_cluster_matrix = clinical_features_matrix[ clin_features , names(which(cluster_assignments == cluster_id)) ]
    }
    if (COUNTER !=0){
      patient_cluster_matrix = cbind(patient_cluster_matrix, clinical_features_matrix[ clin_features , names(which(cluster_assignments == cluster_id)) ] )
    }
    COUNTER = COUNTER + 1
  }
  
  #arbitrary clusters lets you create your own clusters instead of using the clusters from a cluster object  
  if (is.null(arbitrary_clusters) == FALSE){
    nb_clusters = length(arbitrary_clusters)
    COUNTER = 0
    for (i in 1:nb_clusters){
      assign(x = paste0('cluster_',i), value = (COUNTER + 1) : (COUNTER + (arbitrary_clusters[i])) )
      print(get(x = paste0('cluster_',i)))
      COUNTER = COUNTER + (arbitrary_clusters[i]) - 1
    }
  }
  
  #If the value "arbitrary_clusters" is NULL (default) we use the clusters assigned in the cluster object
  #Creates variables cluster_1, cluster_2, etc each with "coordinates" for example cluster one could span patients 1:30, then 31:45 for cluster_2, etc
  if (is.null(arbitrary_clusters) == TRUE){
    COUNTER = 0
    for (i in 1:nb_clusters){
      assign(x = paste0('cluster_',i), value = (COUNTER + 1) : (COUNTER + (table(cluster_assignments)[i])) )
      print(get(x = paste0('cluster_',i)))
      COUNTER = COUNTER + (table(cluster_assignments)[i]) - 1
    } 
  }
  
  #Create data frames to report the clinical data from each cluster (the sums, the porportions, the null ratios)  
  df_report_ratios = data.frame()
  df_report_sums = data.frame()
  df_report_null_ratios = data.frame()
  
  for (i in 1:nb_features){
    print(i)
    for (j in 1:nb_clusters){
      print(paste0(row.names(patient_cluster_matrix)[i], ' | Cluster ', j , ' ',  round(sum(as.numeric(as.factor(patient_cluster_matrix[i,  get(x = paste0('cluster_',j))  ])) -1) / sum(as.numeric(as.factor(patient_cluster_matrix[i,  unlist(sapply(1:nb_clusters, function(k) get(x = paste0('cluster_', k))))    ]))-1) ,2) ))
      df_report_sums[i,j] = sum(as.numeric(as.factor(patient_cluster_matrix[i,  get(x = paste0('cluster_',j))  ])) -1)
      df_report_ratios[i,j] = round( sum(as.numeric(as.factor(patient_cluster_matrix[i,  get(x = paste0('cluster_',j))  ])) -1) / sum(as.numeric(as.factor(patient_cluster_matrix[i,  unlist(sapply(1:nb_clusters, function(k) get(x = paste0('cluster_', k))))     ])) -1) , 2)
      df_report_null_ratios[i,j] = round(length(get(x = paste0('cluster_',j))) / length(unlist(sapply(1:nb_clusters, function(i) get(x = paste0('cluster_', i))))) , 2) 
    }
    df_report_sums[i,(nb_clusters + 1)] = sum(as.numeric(as.factor(patient_cluster_matrix[i,    unlist(sapply(1:nb_clusters, function(i) get(x = paste0('cluster_', i))))   ])) -1)
    df_report_ratios[i,(nb_clusters + 1)] = sum(as.numeric(as.factor(patient_cluster_matrix[i,  unlist(sapply(1:nb_clusters, function(i) get(x = paste0('cluster_', i))))   ])) -1)
    df_report_null_ratios[i,(nb_clusters + 1)] = sum(as.numeric(as.factor(patient_cluster_matrix[i,  unlist(sapply(1:nb_clusters, function(i) get(x = paste0('cluster_', i))))   ])) -1)
  }
  
  
  colnames(df_report_ratios) = c(paste0('cluster_',1:nb_clusters), 'n')
  row.names(df_report_ratios) = row.names(patient_cluster_matrix)[1:nb_features]
  colnames(df_report_sums) = c(paste0('cluster_',1:nb_clusters), 'n')
  row.names(df_report_sums) = row.names(patient_cluster_matrix)[1:nb_features]
  
  list_df = list()
  list_df[['ratios']] = df_report_ratios
  list_df[['null_ratios']] = df_report_null_ratios
  list_df[['sums']] = df_report_sums
  
  #Chi-square test of proportions for clinical features falling in the clusters
  #use TryCatch because each row with sum = 0 will interrupt the program with an error (here I make the p-vals = 1 in case of an error)
  last_column_df_report_sums = ncol(df_report_sums)
  for (i in 1:nrow(df_report_sums)){
    
    df_chisq = data.frame(matrix(NA, nrow = 2 , ncol = nb_clusters))
    df_chisq[1,] = df_report_sums[i , 1:nb_clusters]
    
    
    
    df_chisq[2,] = sapply(1:length(table(cluster_assignments)), function(i) table(cluster_assignments)[i] - df_chisq[1,i] )
    
    print(paste0('Contingency matrix for chi-square. Row 1: sum of 1s ; Row 2: sum of 0s :', row.names(df_report_sums)[i]))
    print(df_chisq)
    
    df_report_sums[i,(last_column_df_report_sums + 1)] = tryCatch(expr = chisq.test(df_chisq)$p.val, error=function(e) 1)
  }
  colnames(df_report_sums)[last_column_df_report_sums + 1] = 'p-val'
  
  
  
  #correct the p-values if FDR == TRUE
  if (FDR == TRUE){
    df_report_sums[,'p-val'] = p.adjust(p = df_report_sums[,'p-val'], method = 'fdr', n = length(df_report_sums[,'p-val']))  
  }
  
  list_df[['sums']] = df_report_sums
  
  color_vector = NULL
  color_vector[1:nrow(df_report_sums)] = 'white'
  color_vector[which(df_report_sums[,'p-val'] < p_val_threshold)] = 'green'
  
  #plotting the features
  par(mfrow = c(nb_features,1), mai = c(0,0,0,0))
  
  for (i in 1:nb_features){
    image(x = seq(from = 0.5, to = (length(cluster_assignments) - 0.5), by = 1),
          y = 1,
          z = as.matrix(as.numeric(as.factor(patient_cluster_matrix[i,]))),
          ylim = c(0,0), ylab = '', xlab = '', axes = F, col = c('black', color_vector[i]))
    
    print(paste0(i,': ', row.names(patient_cluster_matrix)[i]))
    
    abline(v = 1:(length(cluster_assignments)), lty = 1, col = 'black')
    abline(h = 1:(nb_features), lty = 1, col = 'white')
    COUNTER = 0
    for (j in 1:nb_clusters){
      COUNTER = COUNTER + table(cluster_assignments)[[j]]
      abline(v = COUNTER, lty = 2, col = 'yellow')
    }
  }
  
  return(list_df)
  
}



LOOCV_viper <- function(eset, eset_ref, network){
  list_stouffers = list()
  list_samples = colnames(eset)
  
  for (i in 1:ncol(eset)){
    print(paste0('Leaving out: ', colnames(eset)[i]))
    print(paste0('Computing ', i, ' of ', ncol(eset), '...'))
    vps = viperSignature(eset[,-i], eset_ref)
    vp = viper(eset = vps$signature, regulon = network, method = 'none')
    stouffer = apply(vp, 1, function(i) sum(i)/sqrt(length(i)))
    stouffer = sort(stouffer, decreasing = T)
    list_stouffers[[i]] = stouffer
  }
  return(list_stouffers)
}



multiple_intersect <- function(...){
  args = list(...)
  res = unlist(args[1])
  for (i in 2:length(args)){
    #print(res); print(args[i])
    res = intersect(res, unlist(args[i]))
  }
  return(unlist(res))
}

raw_counts_transformation <- function(raw_counts_mat, method = c('cpm', 'logcpm', 'double rank')){
  if (method == 'logcpm'){
    res = log2(t(t(raw_counts_mat)/(colSums(raw_counts_mat)/1e6)) + 1)
  }
  return(res)
}


df_rr <- Plot_Clinical_Features(clinical_features_matrix = Clin_features, cluster_assignments = CC_GEM_Diag_pearson[[2]]$consensusClass, p_val_threshold = 0.05, FDR = T)

######################################################################################
#Load the data (RNAseq and Clinical) and run viper

#Data (pediatric AML) original counts data from TARGET 

original_raw_counts_p_ensembl = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/TARGET_AML/Consolidated_datasets/OriginalRawCounts_p_Ensembl.txt'
original_raw_counts_p_ensembl = read.delim(file = original_raw_counts_p_ensembl, header = T, sep = '\t', row.names = 1, as.is = T, stringsAsFactors = F)


substr(colnames(original_raw_counts_p_ensembl), 2, 18)
substr(colnames(vst_p), 8, 24)

sum(sapply(1:358, function(i) substr(colnames(original_raw_counts_p_ensembl), 2, 18)[i] == substr(colnames(vst_p), 8, 24)[i]))


TF <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/ARACNe/TF/TF_1877.txt', header = F, sep = '\n')
coTF <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/ARACNe/TF/CoTF_677.txt', header = F, sep = '\n')
Sig <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/ARACNe/TF/Sig_3739.txt', header = F, sep = '\n')

raw_counts_p <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/TARGET_AML/Consolidated_datasets/TARGET2021full_counts_entrez.txt', header = T, sep = '\t', row.names = 1, as.is = T, stringsAsFactors = F)
plot(hist(colSums(raw_counts_p[,1:303])), col = rgb(0,1,0,0.5), freq = F)
hist(colSums(raw_counts_p[,304:358]), col = rgb(1,0,0,0.5), freq = F, add = T)

plot(hist(colSums(raw_counts_p[,1:100])), col = rgb(0,1,0,0.5), freq = F)
hist(colSums(raw_counts_p[,101:200]), col = rgb(1,0,0,0.5), freq = F, add = T)

plot(hist(colSums(vst_p[,1:303])), col = rgb(0,1,0,0.5), freq = F)
hist(colSums(vst_p[,304:358]), col = rgb(1,0,0,0.5), freq = F, add = T)

vst_p <- '/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/R/AML_Phenotype_Prediction/Input/GEM_pAML_vst.txt'
vst_p <- read.delim(file = vst_p, header = T, sep = '\t', row.names = 1)
vst_p <- as.matrix(vst_p)

regulon_path <- '/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/R/AML_Phenotype_Prediction/Input/pAML_Network.txt' #Path to the regulon computed by ARACNe
pAML <- aracne2regulon(afile = regulon_path, eset = vst_p, format = '3col', verbose = T)
pAML_TF_coTF <- '/Users/Alexandre/Desktop/Columbia/CalifanoLab/ARACNe/ReconstructedNetworks/TARGET_AML/Ped_AML_2021_Normalized_TF_coTF.txt'
pAML_TF_coTF <- aracne2regulon(afile = pAML_TF_coTF, eset = vst_p, format = '3col', verbose = T)




Clinical <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/R/AML_Phenotype_Prediction/Input/Clinical_Data.txt', header = T, sep = '\t', stringsAsFactors = F)
Clinical$Cohort <- 'NA'
Clinical$Cohort[grep(pattern = '20', x = Clinical$TARGET.USI)] <- 'TARGET.20'
Clinical$Cohort[grep(pattern = '21', x = Clinical$TARGET.USI)] <- 'TARGET.21'
Clinical$TARGET.USI <- gsub(pattern = '.*-', replacement = '', x = Clinical$TARGET.USI)

Clinical$Patient <- substr(Clinical$TARGET.USI, 11, 16)


#Create vectors of patient names, relapse, induction failure and non relapse
#criteria for non-relapse: censored + event-free survival more than 3 years
Patients_T20_diag <- colnames(vp_T20_diag)
Patients_T20_diag <- gsub(pattern = 'TARGET.20.', replacement = '', x = Patients_T20_diag)
Patients_T20_diag <- gsub(pattern = '\\..*.', replacement = '', x = Patients_T20_diag)
relapse <- Clinical$TARGET.USI[Clinical$First.Event == 'Relapse' & Clinical$Event.Free.Survival.Time.in.Days < 365 & Clinical$TARGET.USI %in% Patients_T20_diag]
non.relapse <- Clinical$TARGET.USI[Clinical$First.Event == 'Censored' & Clinical$Event.Free.Survival.Time.in.Days > 1095 & Clinical$TARGET.USI %in% Patients_T20_diag]
relapse <- colnames(vp_T20_diag)[which(Patients_T20_diag %in% relapse)]
non.relapse <- colnames(vp_T20_diag)[which(Patients_T20_diag %in% non.relapse)]
first.events <- sapply(1:ncol(vst_p), function(i) print(Clinical$First.Event[Clinical$TARGET.USI ==  substr(colnames(vst_p),11,16)[i]  ]) )
induction.failure <- colnames(vst_p)[which(first.events == 'Induction failure')]
induction.failure <- induction.failure[grep(pattern = paste0(c('09A', '03A'), collapse = '|'), x = induction.failure )]


#Select patients by FAB subtype
FAB_NR <- sapply(  X = substr(x = non.relapse, 11, 16), FUN = function(i) print(Clinical$FAB.Category[Clinical$TARGET.USI == i]  ))
FAB_R <- sapply(  X = substr(x = relapse, 11, 16), FUN = function(i) print(Clinical$FAB.Category[Clinical$TARGET.USI == i]  ))
FAB_NR_R <- sapply(  X = substr(x = c(non.relapse, relapse), 11, 16), FUN = function(i) print(Clinical$FAB.Category[Clinical$TARGET.USI == i]  ))
FAB_NR_R[FAB_NR_R == 'M0']
c(non.relapse,relapse)[FAB_NR_R == 'M1']




#Identification of matched samples



#VIPER matrices
#vp_full (T20 & T21), vp_T20, vp_T20_diag
vp_p <- viper(eset = vst_p, regulon = pAML_TF_coTF, method = 'ttest', verbose = T)

Diagnostic_T20 <- grep(pattern = paste(c('09A', '03A'),collapse = '|'), x = colnames(vst_p))
Diagnostic_T20 <- colnames(vst_p)[Diagnostic_T20]
Diagnostic_T20 <- Diagnostic_T20[grep(pattern = '20', x = Diagnostic_T20)]

vp_T20_diag <- viper(eset = vst_p[,Diagnostic_T20], regulon = pAML, method = 'ttest', verbose = T)
vp_T20_diag_TF_coTF <- viper(eset = vst_p[,Diagnostic_T20], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)

vp_T20_NR_R_TF_coTF <- viper(eset = vst_p[,c(non.relapse, relapse)], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)


vp_T20T21_Induction_Diag <- viper(eset = vst_p[,c(non.relapse, induction.failure)], regulon = pAML, method = 'ttest', verbose = T)
vp_T20T21_Induction_Diag_TF_coTF <- viper(eset = vst_p[,c(non.relapse, induction.failure)], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)

###############################################################################################################
#VIPER matrices for FAB subtypes M1, M2, M4, M5
####################---------------M2

vp_M1 <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M1']], regulon = pAML, method = 'ttest', verbose = T)

# vps_M2 <- viperSignature(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M2']], 
#                          as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  ,  c(non.relapse,relapse)[FAB_NR_R == 'M2'])], 
#                          per = 2500)

# vp_M2 <- viper(eset = vps_M2, regulon = pAML, method = 'ttest', verbose = T)

vp_M2 <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M2']], regulon = pAML, method = 'ttest', verbose = T)

####################---------------M4
# vps_M4 <- viperSignature(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M4']], 
#                          as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  ,  c(non.relapse,relapse)[FAB_NR_R == 'M4'])], 
#                          per = 2500)

# vp_M4 <- viper(eset = vps_M4, regulon = pAML, method = 'ttest', verbose = T)

vp_M4 <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M4']], regulon = pAML, method = 'ttest', verbose = T)

####################---------------M5
# vps_M5 <- viperSignature(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M5']], 
#                          as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  ,  c(non.relapse,relapse)[FAB_NR_R == 'M5'])], 
#                          per = 2500)

# vp_M5 <- viper(eset = vps_M5, regulon = pAML, method = 'ttest', verbose = T)

vp_M5 <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M5']], regulon = pAML, method = 'ttest', verbose = T)


####################---------------M0_M1_M6_M7_NOS_Unknown
# vps_M0_M1_M6_M7_NOS_Unknown <- viperSignature(eset = as.matrix(vst_p)[, c(non.relapse,relapse) [FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'] ], 
#                                               as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  , c(non.relapse,relapse)[FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'])], 
#                                               per = 2500)

# vp_M0_M1_M6_M7_NOS_Unknown <- viper(eset = vps_M0_M1_M6_M7_NOS_Unknown, regulon = pAML, method = 'ttest', verbose = T)

vp_M0_M1_M6_M7_NOS_Unknown <- viper(eset = as.matrix(vst_p)[, c(non.relapse,relapse) [FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'] ], 
                                                                                   as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  , c(non.relapse,relapse)[FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'])], 
                                    regulon = pAML, method = 'ttest', verbose = T)



#VIPER matrices for subtypes only with TFs and coTFs (disregarding signaling molecules)
vp_M1_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M1']], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)
vp_M2_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M2']], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)
vp_M4_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M4']], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)
vp_M5_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'M5']], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)

vp_M0_M1_M6_M7_NOS_Unknown_TF_co_TF <- viper(eset = as.matrix(vst_p)[, c(non.relapse,relapse) [FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'] ], 
                                    as.matrix(vst_p)[,setdiff ( c(non.relapse,relapse)  , c(non.relapse,relapse)[FAB_NR_R == 'M0' | FAB_NR_R == 'M1' | FAB_NR_R == 'M6' | FAB_NR_R == 'M7' | FAB_NR_R == 'NOS' | FAB_NR_R == 'Unknown'])], 
                                    regulon = pAML_TF_coTF, method = 'ttest', verbose = T)


vp_NOS_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[FAB_NR_R == 'NOS']], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)

vp_M0_M1_M2_TF_coTF <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[c(which(FAB_NR_R == 'M0'), which(FAB_NR_R == 'M1'), which(FAB_NR_R == 'M2'))]], regulon = pAML_TF_coTF, method = 'ttest', verbose = T)
#reorder the columns to have non-relapse as the first samples and relapse as the remaining samples (as opposed to M1-NR-R then M2-NR-R)
vp_M0_M1_M2_TF_coTF <- vp_M0_M1_M2_TF_coTF[,c(which(substr(colnames(vp_M0_M1_M2_TF_coTF), 11, 16) %in% names(FAB_NR)), 
                                              which(substr(colnames(vp_M0_M1_M2_TF_coTF), 11, 16) %in% names(FAB_R)))] #samples 1:39 Non-Relapse; 40:63 Relapse

response_list_M0_M1_M2 <- c(rep(0,39), rep(1,24))

vp_M0_M1_M2 <- viper(eset = as.matrix(vst_p)[,c(non.relapse,relapse)[c(which(FAB_NR_R == 'M0'), which(FAB_NR_R == 'M1'), which(FAB_NR_R == 'M2'))]], regulon = pAML, method = 'ttest', verbose = T)
vp_M0_M1_M2 <- vp_M0_M1_M2[,c(which(substr(colnames(vp_M0_M1_M2), 11, 16) %in% names(FAB_NR)), 
                                              which(substr(colnames(vp_M0_M1_M2), 11, 16) %in% names(FAB_R)))] 

############-----------ADULT AML------------------###############
#ADULT AML
vst_a <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/TCGA/Adult_TCGA_Entrez_Normalized_for_ARACNe.txt', header = T, sep = '\t', row.names = 1, stringsAsFactors = F)
vps_a <- viperSignature(vst_a, vst_a, per = 1000)
vp_a <- viper(eset = vst_a, regulon = aAML, method = 'ttest', verbose = T)

Clinical_a <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/ClinicalData_AdultPed_AML.txt', header = T, sep = '\t', stringsAsFactors = F)
Clinical_a$Patient.ID <- gsub(pattern = '-',replacement = '.', x =  Clinical_a$Patient.ID)


#####db-gap RNA Seq data set ###
dataset_mRNA <- read.delim(file = '/Users/Alexandre/Downloads/SraRunTable.txt', header = T, sep = '\t')

#Totally new samples, not in the public data
dataset_mRNA$Sample_Name[-which(dataset_mRNA$Sample_Name %in% colnames(original_raw_counts_p_ensembl))]
dataset_mRNA$Run[-which(dataset_mRNA$Sample_Name %in% colnames(original_raw_counts_p_ensembl))]

print(paste0(round(sum(dataset_mRNA$MBytes[-which(dataset_mRNA$Sample_Name %in% colnames(original_raw_counts_p_ensembl))])/1024, 1), ' GB'))

New_RAW_counts <- read.delim(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/Raw_data_AML/Final_STAR_Results/TARGET_AML_ENSEMBL_STAR_original.txt', header = T, sep = '\t', row.names = 1, as.is = T)
New_RAW_counts_PatientNames <- New_RAW_counts
colnames(New_RAW_counts_PatientNames) <- sapply(1:419, function(i) dataset_mRNA$Sample_Name[ dataset_mRNA$Run == colnames(New_RAW_counts)[i] ] )

TARGET_21 <- New_RAW_counts_PatientNames[,grep(pattern = 'TARGET-21', x = colnames(New_RAW_counts_PatientNames))]
TARGET_20 <- New_RAW_counts_PatientNames[,grep(pattern = 'TARGET-20', x = colnames(New_RAW_counts_PatientNames))]

New_RAW_counts_PatientNames <- cbind(TARGET_20 , TARGET_21)
write.table(file = '/Users/Alexandre/Desktop/Columbia/CalifanoLab/AML_data/Raw_data_AML/Final_STAR_Results/TARGET_AML_ENSEMBL_STAR.txt', x = New_RAW_counts_PatientNames, quote = F, sep = '\t', row.names = T)
