#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Executes Unit-Tests for Experiment 2.

#########################################################
#' Tests if both formulas lead to the same results.
test_that("Test1", {
  # function 1
  f1 <- function(cluster) {
    summation <- 0;
    
    for (vec1 in cluster) {
      for	(vec2 in cluster) {
        dist <- sum((vec1 - vec2) ^ 2);
        summation <- summation + dist;
      }
    }
    
    return(summation);
  }
  
  # function 2
  f2 <- function(cluster) {
    summation <- 0;
    
    n <- length(cluster);
    average <- Reduce("+", cluster)/n;  # average cluster vector
    
    for (vec in cluster) {
      dist <- sum((vec - average) ^ 2);
      summation <- summation + dist;
    }
    
    return(2*n*summation);
  }
  
  cluster1 <- list(c(1,2), c(2,3), c(3,5));
  cluster2 <- list(c(4,2), c(5,9), c(4,5), c(1,2));
  cluster3 <- list(c(3,2,5), c(5,9,3), c(4,5,2));
  
  testthat::expect_equal(f1(cluster1), f2(cluster1));
  testthat::expect_equal(f1(cluster2), f2(cluster2));
  testthat::expect_equal(f1(cluster3), f2(cluster3));
});


