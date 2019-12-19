Knn_metaCells <- function(normData, rawCounts, K, normalization = c('lognormalize', 'cpm', 'none'), test_k = FALSE){
  require(FNN)
  final_time = 0
  
  if (test_k == TRUE) {
    normalization = 'none'
    res = matrix(NA, nrow = length(K), ncol = 2, dimnames = list(K, c('Coverage', 'Number of genes')))
  }
  
  k_values_to_test = K
  counter = 1
  
  while(counter <= length(k_values_to_test)){
    
    if (final_time != 0) {
      total_time = round(final_time - initial_time)
      message(paste0('Time required for k = ', K, ' : ', total_time, ' seconds'))
    }
    
    K = k_values_to_test[counter]
    message(paste0('\nperforming k-nn with k = ', K))
    
    initial_time = as.numeric(Sys.time())
    #Get the KNN using the clustering on normalized expression data
    knn = get.knn(t(normData), K, algorithm = "brute") 
    
    # Create a matrix of the same size of expmat0
    sum_knn = matrix(0,length(normData[,1]),length(normData[1,]))
    mdata_c1 = rawCounts[ , colnames(normData)]
    
    
    ## check that  the same  cells are selected
    pb = txtProgressBar(1, ncol(mdata_c1), 1, style = 3)
    for (i in 1:ncol(mdata_c1) ){
      sum_knn[,i] = apply(mdata_c1[,colnames(mdata_c1[,c(i,knn$nn.index[i,])])],1,sum)
      setTxtProgressBar(pb, i)
    }
    
    dim(sum_knn)
    rownames(sum_knn) = rownames(mdata_c1)
    colnames(sum_knn) = colnames(mdata_c1)
    
    ind <- colSums(sum_knn)>0
    
    
    if (normalization == 'lognormalize') tmp = log2(t(t(sum_knn[, ind])/(colSums(sum_knn[, ind])/1e6)) + 1)
    if (normalization == 'cpm') tmp = t(t(sum_knn[, ind])/(colSums(sum_knn[, ind])/1e6))
    if (normalization == 'none') tmp = sum_knn
    
    if (test_k == TRUE){
      res[counter, 1] = mean(apply(tmp, 2, sum))
      res[counter, 2] = mean(apply(tmp, 2, function(i) sum(i != 0) ))
    }
    else{
      res = tmp
    }
    counter = counter + 1
    final_time = as.numeric(Sys.time())
  }
  
  
  
  return(res)
}
