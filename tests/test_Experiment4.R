#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Executes Unit-Tests for Experiment 4.

#########################################################
#' Tests if it is the correct implementation.
test_that("Test1", {
  value <- (21-5)/sqrt(27*27);
  correalationCoefficient <- cor(c(8,6,5,3.5,1,2,3.5,7), 
                                 c(6,7.5,4,1,2,3,5,7.5), 
                                 method = "kendall");
  
  testthat::expect_equal(round(value, 4), 
                         round(correalationCoefficient, 4));
});


#########################################################
#' Tests if Tukey's Biweight Robust Mean Function is implemented how desired.
test_that("Test2", {
  result <- 3.288936942;  # by hand computed result (see pdf)
  resultImplementation <- Math.tukeysBiweightRobustMean(c(2,3,5));
  
  testthat::expect_equal(round(result, 4), 
                         round(resultImplementation, 4));
});