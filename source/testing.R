testing <- function(mean_matrix, sd_matrix, hypo_meth_ranks, hyper_meth_ranks){
  mean_test <- mean_matrix[12345,]
  sd_test <- sd_matrix[12345,]
  mean_1_sd_test <- mean_test - sd_test
  mean_2_sd_test <- mean_test - 2 * sd_test
  
  #test hypo methylation results
  testthat::expect_equal(hypo_meth_ranks$`number of cell types with score lower than average`[12345,],
                         (rank(mean_test, ties.method = "max" )-1))
  
  testthat::expect_equal(hypo_meth_ranks$`number of cell types with score lower than (average - SD)`[12345,5],
                         (rank(c(mean_test[5], mean_1_sd_test[-5]), ties.method = "max" )-1)[1])
  
  testthat::expect_equal(hypo_meth_ranks$`number of cell types with score lower than (average - 2*SD)`[12345,35],
                         (rank(c(mean_test[35], mean_2_sd_test[-35]), ties.method = "max" )-1)[1])
  
  #test hyper methylation results
  testthat::expect_equal(hyper_meth_ranks$`number of cell types with score higher than average`[12345,],
                         (rank(-mean_test, ties.method = "max" )-1))
  
  testthat::expect_equal(hyper_meth_ranks$`number of cell types with score higher than (average + SD)`[12345,5],
                         (rank(c(-mean_test[5], -mean_1_sd_test[-5]), ties.method = "max" )-1)[1])
  
  testthat::expect_equal(hyper_meth_ranks$`number of cell types with score higher than (average + 2*SD)`[12345,35],
                         (rank(c(-mean_test[35], -mean_2_sd_test[-35]), ties.method = "max" )-1)[1])
}