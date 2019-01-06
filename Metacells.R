Knn_metaCells_Pasquale <- function(normData, rawCounts, K)
{
  require(FNN)
  #Get the KNN using the clustering onnormalized expression data
  knn_10<-get.knn(t(normData),K,algorithm = "brute") 
  
  # Create a matrix of the same size of expmat0
  sum_knn<-matrix(0,length(normData[,1]),length(normData[1,]))
  mdata_c1<-rawCounts[,colnames(normData)]
  
  message("Check size")
  dim(mdata_c1)
  
  ## check that  the same  cells are selected
  for (i in 1:length(mdata_c1[1,]))
  {
    sum_knn[,i]<-apply(mdata_c1[,colnames(mdata_c1[,c(i,knn_10$nn.index[i,])])],1,sum)
    print(paste0("MetaCell_",i))
  }
  dim(sum_knn)
  rownames(sum_knn)<-rownames(mdata_c1)
  colnames(sum_knn)<-colnames(mdata_c1)
  
  ind <- colSums(sum_knn)>0
  ### Now let's compute tpm, signature and protein activity
  tpm_KNN <- log2(t(t(sum_knn[, ind])/(colSums(sum_knn[, ind])/1e6)) + 1)
  return(tpm_KNN)
}

metacell_pasquale_diag <- Knn_metaCells_Pasquale(normalized_sc_data[,grep('\\.1', colnames(normalized_sc_data))], sparse_raw_counts[,grep('\\.1', colnames(sparse_raw_counts))], 10)
metacell_pasquale_relapse <- Knn_metaCells_Pasquale(normalized_sc_data[,grep('\\.2', colnames(normalized_sc_data))], sparse_raw_counts[,grep('\\.2', colnames(sparse_raw_counts))], 10)
#metacell_all_cells <- Knn_metaCells(metacell_all_cells, raw_counts_all_cells, 10)

cpm_metacells_pasquale_diag <- raw_counts_transformation(metacell_pasquale_diag, 'logcpm')
cpm_metacells_pasquale_relapse <- raw_counts_transformation(metacell_pasquale_relapse, 'logcpm')

PCA(cpm_metacells_pasquale_diag[order(rowVars(cpm_metacells_pasquale_diag), decreasing = T)[1:500] , sample(colnames(cpm_metacells_pasquale_diag), round(0.1*ncol(cpm_metacells_pasquale_diag)), F)])
PCA(cpm_metacells_pasquale_relapse[order(rowVars(cpm_metacells_pasquale_relapse), decreasing = T)[1:500] , sample(colnames(cpm_metacells_pasquale_relapse), round(1*ncol(cpm_metacells_pasquale_relapse)), F)])


write.table(x = metacell_pasquale_diag, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/ALL/FerrandoData/scRNAseq/scRNAseq_DxRe_TALL/metacell_pasquale_diag.txt', quote = F, sep = '\t', row.names = T)
write.table(x = metacell_pasquale_relapse, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/ALL/FerrandoData/scRNAseq/scRNAseq_DxRe_TALL/metacell_pasquale_relapse.txt', quote = F, sep = '\t', row.names = T)
#write.table(x = metacell_all_cells, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/ALL/FerrandoData/scRNAseq/scRNAseq_DxRe_TALL/metacell_all_cells.txt', quote = F, sep = '\t', row.names = T)

library(Seurat)
seurat_object_sc_T_ALL <- Read10X('/Users/Alexandre/Desktop/Columbia/CalifanoLab/data/ALL/FerrandoData/scRNAseq/scRNAseq_DxRe_TALL/filtered_gene_bc_matrices/GRCh38/')
seurat_object_sc_T_ALL <- CreateSeuratObject(seurat_object_sc_T_ALL, min.cells = 3, min.genes = 200, project = 'ALL_one_patient')

seurat_object_sc_T_ALL@raw.data <- sparse_raw_counts
seurat_object_sc_T_ALL@data <- normalized_sc_data
seurat_object_sc_T_ALL <- FindVariableGenes(object = seurat_object_sc_T_ALL, mean.function = ExpMean, dispersion.function = LogVMR, 
                                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

mito.genes <- grep(pattern = "^MT-", x = row.names(x = seurat_object_sc_T_ALL@data), value = TRUE)
percent.mito <- colSums(seurat_object_sc_T_ALL@raw.data[mito.genes, ]) / colSums(seurat_object_sc_T_ALL@raw.data)

colnames(sparse_raw_counts) <- gsub('\\.', '-', colnames(sparse_raw_counts))
colnames(normalized_sc_data) <- gsub('\\.', '-', colnames(normalized_sc_data))

colnames(seurat_object_sc_T_ALL@raw.data) <- gsub('\\.', '-', colnames(seurat_object_sc_T_ALL@raw.data))
colnames(seurat_object_sc_T_ALL@data) <- gsub('\\.', '-', colnames(seurat_object_sc_T_ALL@data))

seurat_object_sc_T_ALL@ident <- seurat_object_sc_T_ALL@ident[colnames(seurat_object_sc_T_ALL@raw.data)]
seurat_object_sc_T_ALL@meta.data <- seurat_object_sc_T_ALL@meta.data[colnames(seurat_object_sc_T_ALL@raw.data),]

seurat_object_sc_T_ALL@spatial@mix.probs[colnames(seurat_object_sc_T_ALL@raw.data)]

seurat_object_sc_T_ALL@meta.data$percent.mito <- percent.mito
seurat_object_sc_T_ALL@cell.names <- colnames(seurat_object_sc_T_ALL@raw.data)
seurat_object_sc_T_ALL <- ScaleData(object = seurat_object_sc_T_ALL, vars.to.regress = c("nUMI", "percent.mito"))
seurat_object_sc_T_ALL <- RunPCA(object = seurat_object_sc_T_ALL, pc.genes = seurat_object_sc_T_ALL@var.genes, do.print = TRUE, pcs.print = 1:5, 
                                 genes.print = 5)

PCAPlot(object = seurat_object_sc_T_ALL, dim.1 = 1, dim.2 = 2)
seurat_object_sc_T_ALL <- FindClusters(object = seurat_object_sc_T_ALL, reduction.type = "pca", save.SNN = TRUE)

metacell_ajay_diag <- Knn_metaCells_Ajay(seurat_object = seurat_object_sc_T_ALL, 
                                         normalized_counts = normalized_sc_data[,grep('-1', colnames(normalized_sc_data))], 
                                         raw_counts = sparse_raw_counts[,grep('-1', colnames(sparse_raw_counts))], 
                                         k = 10, n = 2)

PCA(metacell_ajay_diag$metacell_raw_counts)
PCA(raw_counts_transformation(metacell_ajay_diag$metacell_raw_counts, 'logcpm'), show_legend = F)

metacell_ajay_relapse <- Knn_metaCells_Ajay(seurat_object = seurat_object_sc_T_ALL, 
                                         normalized_counts = normalized_sc_data[,grep('-2', colnames(normalized_sc_data))], 
                                         raw_counts = sparse_raw_counts[,grep('-2', colnames(sparse_raw_counts))], 
                                         k = 10, n = 2)

PCA(raw_counts_transformation(metacell_ajay_relapse$metacell_raw_counts, 'logcpm'), show_legend = F)

merged_metacells_diag_relapse_ajay <- cbind(metacell_ajay_diag$metacell_raw_counts, metacell_ajay_relapse$metacell_raw_counts)
PCA(raw_counts_transformation(merged_metacells_diag_relapse_ajay, 'logcpm'), show_legend = F, pt_colors = c(rep('blue', 166),
                                                                                                            rep('red', 182)))




metacell_ajay_diag_relapse <- Knn_metaCells_Ajay(seurat_object = seurat_object_sc_T_ALL, 
                                                 normalized_counts = normalized_sc_data[,grep('-2', colnames(normalized_sc_data))], 
                                                 raw_counts = sparse_raw_counts[,grep('-2', colnames(sparse_raw_counts))], 
                                                 k = 10, n = 2)

#diff gene experssion between 
#number of genes detected vs k



# entrez_metacell_cluster_1 <- convert_genenames_matrix(metacell_cluster_1, 'genenames', 'entrez')
# entrez_metacell_cluster_2 <- convert_genenames_matrix(metacell_cluster_2, 'genenames', 'entrez')
# entrez_all_cells <- convert_genenames_matrix(metacell_all_cells, 'genenames', 'entrez')

TF <- as.character(read.delim('/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/TF.txt', header = F, sep = '\t', as.is = T)[[1]])
coTF <- as.character(read.delim('/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/coTF.txt', header = F, sep = '\t', as.is = T)[[1]])
sig <- as.character(read.delim('/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/sig.txt', header = F, sep = '\t', as.is = T)[[1]])
surf <- as.character(read.delim('/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/surface_markers.txt', header = F, sep = '\t', as.is = T)[[1]])

TF_genenames <- Genename_Conv('entrez', 'genenames', TF) ; TF_genenames <- TF_genenames[-which(is.na(TF_genenames))]
coTF_genenames <- Genename_Conv('entrez', 'genenames', coTF) ; coTF_genenames <- coTF_genenames[-which(is.na(coTF_genenames))]
sig_genenames <- Genename_Conv('entrez', 'genenames', sig) ; sig_genenames <- sig_genenames[-which(is.na(sig_genenames))]
surf_genenames <- Genename_Conv('entrez', 'genenames', surf) ; surf_genenames <- surf_genenames[-which(is.na(surf_genenames))]

write.table(x = TF_genenames, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/TF_genenames.txt', quote = F, sep = '\n', row.names = F, col.names = F)
write.table(x = coTF_genenames, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/coTF_genenames.txt', quote = F, sep = '\n', row.names = F, col.names = F)
write.table(x = sig_genenames, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/sig_genenames.txt', quote = F, sep = '\n', row.names = F, col.names = F)
write.table(x = surf_genenames, file = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/aracne/lists/surf_genenames.txt', quote = F, sep = '\n', row.names = F, col.names = F)


sum(row.names(entrez_all_cells) %in% TF)
sum(row.names(entrez_all_cells) %in% coTF)
sum(row.names(entrez_all_cells) %in% sig)
sum(row.names(entrez_all_cells) %in% surf)

sum(row.names(metacell_all_cells) %in% TF_genenames)
sum(row.names(metacell_all_cells) %in% coTF_genenames)
sum(row.names(metacell_all_cells) %in% sig_genenames)
sum(row.names(metacell_all_cells) %in% surf_genenames)

