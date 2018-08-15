#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Executes Unit-Tests for Experiment 3.

#########################################################
#' Tests if simple weighting works.
test_that("Test1", {
  bestYearsPerColumn <- list();
  
  weighting <- 3;
  bestYears <- c(1992, 1502, 1493);
  
  bestYearsPerColumn <- Analyzer.__getWeightedYearsPerColumn(bestYearsPerColumn, weighting, 
                                                             bestYears, FALSE);
  testthat::expect_equal(bestYearsPerColumn, 
                         list(c(1992, 1502, 1493), c(1992, 1502, 1493), c(1992, 1502, 1493)));
});


#########################################################
#' Tests if the double weighting works.
test_that("Test2", {
  bestYearsPerColumn <- list();
  
  weighting <- 3;
  bestYears <- c(1992, 1502, 1493, 1801, 1723);

  bestYearsPerColumn <- Analyzer.__getWeightedYearsPerColumn(bestYearsPerColumn, weighting, 
                                                             bestYears, TRUE);

  testthat::expect_equal(bestYearsPerColumn, 
                         list(1992, 1992, 1992, 1502, 1502, 1493, 1801, 1723));
});


#########################################################
#' Tests if peak selection works.
test_that("Test3", {
  t1 <- list(list(1992, 1502, 1493), list(1992, 1502, 1502));
  t2 <- list(list(1945, 1484, 1954), list(1945, 1484, 1954));
  t3 <- list(list(2020, 1253, 1832), list(2020, 1832, 1832));
  topYearsPerSample <- list(t1, t2, t3);
  
  samplesYears1 <- list(years = 1992:1997);
  sampleYears2 <- list(years = 1484:1489);
  sampleYears3 <- list(years = 1832:1837);
  
  sampleData <- list(samplesYears1, sampleYears2, sampleYears3);

  ranks <- Analyzer.getPeakRanks(sampleData, topYearsPerSample);

  testthat::expect_equal(ranks, list(2, 3, 1));
});


#########################################################
#' Tests if power set table creation works.
test_that("Test4", {
  scoresTables <- Loader.readInIndividualScores(Paths.TEST_SCORES_2, 
                                               Symbols.REAL_DATA_SEPARATOR)$scoresTables;
  
  powersetTableData <- Analyzer.__getPowerSetTable(scoresTables[[1]]);
  scoresTable <- powersetTableData$scoresTable;
  combinations <- powersetTableData$combinations;
  
  years <- scoresTable[,1];
  scores1 <- scoresTable[,2];
  scores2 <- scoresTable[,3];
  scores3 <- scoresTable[,4];
  scores4 <- scoresTable[,5];
  scores5 <- scoresTable[,6];
  scores6 <- scoresTable[,7];
  scores7 <- scoresTable[,8];
  
  testthat::expect_equal(years, c(1992, 1993, 1994));
  testthat::expect_equal(scores1, c(0.1, 0.4, 0.7));
  testthat::expect_equal(scores2, c(0.2, 0.5, 0.8));
  testthat::expect_equal(scores3, c(0.3, 0.6, 0.9));
  testthat::expect_equal(scores4, c(0.3, 0.9, 1.5));
  testthat::expect_equal(scores5, c(0.4, 1.0, 1.6));
  testthat::expect_equal(scores6, c(0.5, 1.1, 1.7));
  testthat::expect_equal(scores7, c(0.6, 1.5, 2.4));
  
  testthat::expect_equal(length(combinations), 7);
});
