#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Executes Unit-Tests for Experiment 1.

#########################################################
#' Tests if the conversion to pure coordinates works.
test_that("Test1", {
  absoluteValues <- list(5, 3, 2);
  pureValues <- Math.computePureValues(absoluteValues);

  testthat::expect_equal(round(1.09109, digits=2), round(pureValues[[1]], digits=2));
  testthat::expect_equal(round(-0.2182179, digits=2), round(pureValues[[2]], digits=2));
  testthat::expect_equal(round(-0.8728717, digits=2), round(pureValues[[3]], digits=2));
});


#########################################################
#' Tests if the mean distance computation works.
test_that("Test2", {
  a <- list(2, 3, 5);
  b <- list(3, 1, 2);
  distance <- Math.getMeanDistance(a, b);
  
  testthat::expect_equal(2, distance);
});


#########################################################
#' Test ROC-Plots.
test_that("Test3", {
  ranks <- list(1, 1, 1, 2, 3, 5, 6, 1, 1, 4);
  scores <- list(0.04, 0.08, 0.11, 0.08, 0.22, 0.31, 0.12, 0.15, 0.04, 0.2);
  ranksPerScoreType <- list(ranks);
  scoresPerScoreType <- list(scores);
  ranksPerPass <- list(ranksPerScoreType);
  scoresPerPass <- list(scoresPerScoreType);

  passes <- list(1);
  scoreTypes <- list("test");

  combinedDataPerScoreType <- Encoder.createCombinedData(ranksPerPass, scoresPerPass, list(1), scoreTypes);

  rates <- Plotter.plotROC(combinedDataPerScoreType,
                           scoreTypes, Titles.CHARACTERISTIC_3,
                           TRUE, FALSE, Extensions.PNG, 3);

  testthat::expect_equal(0, rates$truePositives[[1]]);  # d = 0.04
  testthat::expect_equal(0, rates$falsePositives[[1]]);

  testthat::expect_equal(0, rates$truePositives[[1]]);  # d = 0.04
  testthat::expect_equal(0, rates$falsePositives[[1]]);

  testthat::expect_equal(round(2/7, 3), round(rates$truePositives[[3]], 3));  # d = 0.08
  testthat::expect_equal(0, rates$falsePositives[[3]]);

  testthat::expect_equal(round(2/7, 3), round(rates$truePositives[[4]], 3));  # d = 0.08
  testthat::expect_equal(0, rates$falsePositives[[4]]);

  testthat::expect_equal(round(4/7, 3), round(rates$truePositives[[5]], 3));  # d = 0.11
  testthat::expect_equal(0, rates$falsePositives[[5]]);

  testthat::expect_equal(round(5/7, 3), round(rates$truePositives[[6]], 3));  # d = 0.12
  testthat::expect_equal(0, rates$falsePositives[[6]]);

  testthat::expect_equal(round(5/7, 3), round(rates$truePositives[[7]], 3));  # d = 0.15
  testthat::expect_equal(round(1/3, 3), round(rates$falsePositives[[7]], 3));

  testthat::expect_equal(round(6/7, 3), round(rates$truePositives[[8]], 3));  # d = 0.2
  testthat::expect_equal(round(1/3, 3), round(rates$falsePositives[[8]], 3));

  testthat::expect_equal(round(6/7, 3), round(rates$truePositives[[9]], 3));  # d = 0.22
  testthat::expect_equal(round(2/3, 3), round(rates$falsePositives[[9]], 3));

  testthat::expect_equal(1, round(rates$truePositives[[10]], 3));  # d = 0.31
  testthat::expect_equal(round(2/3, 3), round(rates$falsePositives[[10]], 3));
});


#########################################################
#' Test violin plots.
test_that("Test4", {
  ranks <- list(1, 1, 1, 2, 3, 5, 6, 1, 1, 4);
  ranksPerScoreType <- list(ranks);
  ranksPerPass <- list(ranksPerScoreType);

  passes <- list(1); 
  scoreTypes <- list("test");
  
  combinedRanksPerScoreType <- Encoder.createCombinedRanks(ranksPerPass, passes, scoreTypes);
  
  properties <- Plotter.plotViolinePlots(combinedRanksPerScoreType, scoreTypes, Titles.VIOLINE_RANKS_PLOT, 
                                         TRUE, FALSE, Extensions.PDF, Palettes.NATURE);
  
  testthat::expect_equal(1.5, properties$median);
  testthat::expect_equal(2.5, properties$mean);
  testthat::expect_equal(3.61, properties$variance);
});


#########################################################
#' Testing the ranking function for violin and rank plots.
test_that("Test5", {
  ranks <- Analyzer.getRanks(list(0.02, 0.02, 0.04, 0.03), 
                             list(list(0.05, 0.03, 0.04, 0.02, 0.02, 0.12),
                                  list(0.15, 0.02, 0.03, 0.05, 0.04, 0.08),
                                  list(0.07, 0.18, 0.04, 0.02, 0.05, 0.12),
                                  list(0.06, 0.03, 0.04, 0.17, 0.02, 0.12)));
  
  testthat::expect_equal(2, ranks[[1]]);
  testthat::expect_equal(1, ranks[[2]]);
  testthat::expect_equal(2, ranks[[3]]);
  testthat::expect_equal(2, ranks[[4]]);
});


#########################################################
#' Testing the scores loading function for violin and rank plots.
test_that("Test6", {
  sampleData <- Loader.readInCurves(Paths.TEST_CURVES, Symbols.REAL_DATA_SEPARATOR);
  scoresPerSample <- Loader.readInScores(Paths.TEST_SCORES, Symbols.REAL_DATA_SEPARATOR)$scoresData; 
  scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData, 1960:1974);

  testthat::expect_equal(0.96, scoresForCorrectPositions[[1]]);
  testthat::expect_equal(1.05, scoresForCorrectPositions[[2]]);
  testthat::expect_equal(0.82, scoresForCorrectPositions[[3]]);
});