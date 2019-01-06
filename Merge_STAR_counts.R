merge_RNA_seq_counts_files <- function(PATH, pattern_in_filenames = '\\.txt', stranded = TRUE){
  list_counts = list.files(PATH)
  list_counts = list_counts[grep(pattern = pattern_in_filenames, x = list_counts)]
  
  message('loading the counts files...')
  pb = txtProgressBar(min = 1, max = length(list_counts), initial = 1, style = 3)
  list_res = lapply(list_counts, function(i) {
    setTxtProgressBar(pb, i)
    read.delim(paste0(PATH, '/', i), header = F, sep = '\t') 
  }
    
    )
  names(list_res) = unlist(sapply(list_counts, function(i)  substr(i, 1, nchar(i) - 4))) #remove the ".txt" extension in the sample names
  
  
  
    mat_total_counts = matrix(data = NA, nrow = length(list_res), ncol = 3, dimnames = list(names(list_res), c('unstranded', 'strand_1', 'strand_2')))
    mat_total_counts[1:length(list_res), 1:3] = t(sapply(1:length(list_res), function(i) apply(list_res[[i]][5:nrow(list_res[[1]]), 2:4], 2, sum)))
    mat_total_counts = cbind(mat_total_counts, mat_total_counts[,2]/(mat_total_counts[,2] + mat_total_counts[,3]))
    colnames(mat_total_counts)[4] = 'Ratio'
    mat_total_counts = cbind(mat_total_counts, ifelse(mat_total_counts[,4] < 0.5, 4, 3)) #if ratio of strand 1 to strand 2 is greater than 0.5, column 4 is the column corresponding to the result
    colnames(mat_total_counts)[5] = 'Sequenced_Strand'
    
  
  raw_counts = matrix(data = NA, 
                      nrow = nrow(list_res[[1]]) - 4, 
                      ncol = length(list_res), 
                      dimnames = list(list_res[[1]][5:nrow(list_res[[1]]), 1], names(list_res)  ) 
                      )
  
  message('merging the counts files...')
  
  if (stranded == TRUE){
    raw_counts[1:nrow(raw_counts), 1:ncol(raw_counts)] = sapply(1:length(list_res), 
                                                                function(i) list_res[[i]][5:nrow(list_res[[i]]) , mat_total_counts[i,5]] )  #will pick the strand with highest coverage (column 5 in mat_total_counts)  
  }
  if (stranded == FALSE){
    raw_counts[1:nrow(raw_counts), 1:ncol(raw_counts)] = sapply(1:length(list_res), 
                                                                function(i) {
                                                                  setTxtProgressBar(pb, i)
                                                                  list_res[[i]][5:nrow(list_res[[i]]) , 2]  #will pick the unstranded column (2)
                                                                })
                                                                  
  }
  
 return(raw_counts) 
}





###########
###########
##########
PATH = '/Users/alexandre/Desktop/Columbia/CalifanoLab/data/ALL/dbGaP/raw_counts/T_ALL/consolidated_data/'
list_counts = list.files(PATH)
list_counts = list_counts[grep(pattern = 'txt', x = list_counts)]

list_res = lapply(list_counts, function(i) read.delim(paste0(PATH, '/', i), header = F, sep = '\t') )
names(list_res) = unlist(sapply(list_counts, function(i)  substr(i, 1, nchar(i) - 4)))

list_res[[1]][5:nrow(list_res[[1]]),]
mat_total_counts = matrix(data = NA, nrow = length(list_res), ncol = 3, dimnames = list(names(list_res), c('unstranded', 'strand_1', 'strand_2')))
mat_total_counts[1:length(list_res), 1:3] = t(sapply(1:length(list_res), function(i) apply(list_res[[i]][5:nrow(list_res[[1]]), 2:4], 2, sum)))
mat_total_counts = cbind(mat_total_counts, mat_total_counts[,2]/(mat_total_counts[,2] + mat_total_counts[,3]))
colnames(mat_total_counts)[4] = 'Ratio'
mat_total_counts = cbind(mat_total_counts, ifelse(mat_total_counts[,4] < 0.5, 4, 3)) #if ratio of strand 1 to strand 2 is greater than 0.5, column 4 is the column corresponding to the result
colnames(mat_total_counts)[5] = 'Sequenced_Strand'

View(mat_total_counts)

raw_counts = matrix(data = NA, 
                    nrow = nrow(list_res[[1]]) - 4, 
                    ncol = length(list_res), 
                    dimnames = list(list_res[[1]][5:nrow(list_res[[1]]), 1], names(list_res)  ) 
                    )

raw_counts[1:nrow(raw_counts), 1:ncol(raw_counts)] = sapply(1:length(list_res), 
                                                             function(i) 
                                                               list_res[[i]][5:nrow(list_res[[i]]) , mat_total_counts[i,5]] 
                                                            )

raw_counts[1:5,10]
list_res[[10]][5:10,]

write.table(x = as.matrix(raw_counts), file = paste0(PATH, '/raw_counts_T_ALL_dbGaP.txt'), append = F, quote = F,  sep = '\t', row.names = T, col.names = T)

#merge duplicates
which(duplicated(substr(colnames(raw_counts), 1, 5)))

list_duplicated_samples = unique(sapply(1:ncol(raw_counts), 
                                            function(i) 
                                              which(substr(colnames(raw_counts), 1, 5) == substr(colnames(raw_counts)[i], 1, 5))
                                        )
                                 )

mat_merged_raw_counts = matrix(data = NA, 
                               nrow = nrow(raw_counts), 
                               ncol = length(list_duplicated_samples), 
                               dimnames = list(row.names(raw_counts), unique(substr(colnames(raw_counts), 1, 5))))

for (i in 1:length(list_duplicated_samples)){
  if (length(list_duplicated_samples[[i]]) == 1){
    mat_merged_raw_counts[,i] = raw_counts[, list_duplicated_samples[[i]]]
  }
  
  if (length(list_duplicated_samples[[i]]) == 2){
    mat_merged_raw_counts[,i] = raw_counts[, list_duplicated_samples[[i]][1] ] + raw_counts[, list_duplicated_samples[[i]][2] ]
  }
}

write.table(x = as.matrix(mat_merged_raw_counts), file = paste0(PATH, '/raw_counts_patients_merged_duplicates.txt'), append = F, quote = F,  sep = '\t', row.names = T, col.names = T)


sort(apply(mat_merged_raw_counts, 2, sum))
sort(apply(raw_counts, 2, sum))





AK022_1 <- read.delim(file = paste0(PATH, '/AK022_GAGTGG_L006_001.txt'), header = T, sep = '\t')
AK022_2 <- read.delim(file = paste0(PATH, '/AK022_GAGTGG_L006_002.txt'), header = T, sep = '\t')
apply(AK022_1[4:nrow(AK022_1),2:4], 2, sum)
apply(AK022_2[4:nrow(AK022_2),2:4], 2, sum)





