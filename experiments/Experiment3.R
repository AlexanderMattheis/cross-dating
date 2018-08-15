#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_3_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes third reproducible experiment which was used to create the Voting Approach results.

#########################################################
#' Starts the experiment.
Experiment3.start <- function() {
  # create scores
  # Experiment2.storeComputedScores(Defaults.TASK_A, 3, min, identity);
  # Experiment2.storeComputedScores(Defaults.TASK_B, 3, min, identity);
  # Experiment2.storeComputedScores(Defaults.TASK_C, 3, min, identity);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, identity);
  
  # compute histograms per year to check correctness
  # Experiment3.__computeYearFrequencyHistograms(Paths.PASSES_VOTING, 
  #                                              Files.SCORES_FOLDER_NAME, 1, Defaults.TASK_D, 10);
  
  # compute box plots
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 10, FALSE);
  
  # computes a box plot for the powerset approach
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING_LENGTH_5, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 1, TRUE);
  
  # computes a box plot for the powerset approach with double weighting (i.e. with respect to the ranks)
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING_LENGTH_5, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 4, TRUE, TRUE);
  
  # computes a box plot for the powerset approach with logarithmic weighting (i.e. with respect to the ranks)
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING_LENGTH_5, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 1, TRUE, FALSE, TRUE);
  
  # computes a box plot for the powerset approach and length 10
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 1, TRUE, FALSE, FALSE, 8);
  
  # computes a box plot for the powerset approach and length 15
  # Experiment3.__computeBoxPlots(Paths.PASSES_VOTING_LENGTH_15, Files.SCORES_FOLDER_NAME,
  #                               list(1, 2, 3, 4, 5, 6, 7),  list("a", "b", "c", "d"),
  #                               Palettes.NATURE, 2, TRUE, FALSE, FALSE, 13);
  
  # compute histograms (recomputation took to many time) for length 5 method b)
  # bucket-based, top years 1, top years 2, top years 3
  # Plotter.createRanksHistograms(c(168, 64, 35, 31, 19), c(163, 43, 31, 19, 16),
  #                               c(168, 63, 24, 24, 22), c(155, 68, 40, 23, 10),
  #                               titles = MultiTitles.RANK_HISTOGRAMS, Colors.RED);

  # compute histograms (recomputation took to many time) for length 5 method b)
  # top years 2-5
  # Plotter.createRanksHistograms(c(166, 59, 33, 21, 14), c(167, 58, 37, 21, 12),
  #                               c(163, 62, 41, 20, 15), c(164, 60, 41, 25, 19),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_2, Colors.RED);
  
  # compute histograms (recomputation took to many time) for length 10 methods a-d) min-length 8, top-years 1
  # Plotter.createRanksHistograms(c(97, 21, 9, 7, 7), c(142, 24, 9, 5, 3),
  #                               c(141, 11, 11, 4, 1), c(147, 24, 7, 5, 5),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # compute histograms (recomputation took to many time) for length 15 methods a-d) min-length 8, top-years 1
  # Plotter.createRanksHistograms(c(125, 19, 8, 4, 6), c(161, 19, 4, 7, 4),
  #                               c(166, 9, 4, 3, 4), c(171, 8, 8, 5, 1),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # pure buckets-based approach
  # Plotter.createRanksHistograms(c(126, 20, 12, 4, 7), c(162, 19, 9, 3, 7),
  #                               c(165, 12, 5, 6, 6), c(172, 7, 11, 6, 7),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # evaluate histograms of samples
  # Experiment3.__evaluateHistogramGoodness(passPathGenericName = Paths.PASSES_VOTING,
  #                                         scoresFolderGenericName = Files.SCORES_FOLDER_NAME,
  #                                         passes = list(1, 2, 3, 4, 5),
  #                                         scoreType = Defaults.TASK_D,
  #                                         topYearsCount = 1,
  #                                         minimumLength = 8,
  #                                         xLimits <- c(0, 225),
  #                                         yLimits <- c(-1.1, 1.1));
  
  # plots data for which histograms were used
  # length 5 (t = 1)
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(120, 26, 26, 15, 10), color = Colors.GREEN);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(163, 43, 31, 19, 16), color = Colors.RED);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(159, 41, 20, 16, 15), color = Colors.BLUE);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(151, 48, 40, 14, 11), color = Colors.VIOLET);
  
  # Plotter.createRanksHistograms(c(120, 26, 26, 15, 10), c(163, 43, 31, 19, 16),
  #                               c(159, 41, 20, 16, 15), c(151, 48, 40, 14, 11),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # length 5 (t = 1 up to t = 3)

  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(168, 64, 35, 31, 19), color = Colors.GREEN);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(163, 43, 31, 19, 16), color = Colors.RED);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(168, 63, 24, 24, 22), color = Colors.BLUE);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(155, 68, 40, 23, 10), color = Colors.VIOLET);
  
  # Plotter.createRanksHistograms(c(0, 0, 0, 0, 0), c(163, 43, 31, 19, 16), 
  #                               c(168, 63, 24, 24, 22), c(155, 68, 40, 23, 10), 
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # length 5 double weighted (t = 2 up to t = 4)
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(166, 59, 33, 21, 14), color = Colors.GREEN);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(167, 58, 37, 21, 12), color = Colors.RED);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(163, 62, 41, 20, 15), color = Colors.BLUE);
  
  # Plotter.createRanksHistograms(c(0, 0, 0, 0, 0), c(166, 59, 33, 21, 14), 
  #                               c(167, 58, 37, 21, 12), c(163, 62, 41, 20, 15), 
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # length 10 simple
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(48, 25, 9, 11, 7), color = Colors.GREEN);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(78, 26, 19, 14, 11), color = Colors.RED);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(73, 23, 18, 18, 6), color = Colors.BLUE);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(73, 22, 16, 14, 8), color = Colors.VIOLET);
  
  # Plotter.createRanksHistograms(c(48, 25, 9, 11, 7), c(78, 26, 19, 14, 11),
  #                               c(73, 23, 18, 18, 6), c(73, 22, 16, 14, 8),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  
  # length 10
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(97, 21, 9, 7, 7), color = Colors.GREEN);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(142, 24, 9, 5, 3), color = Colors.RED);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(141, 11, 11, 4, 1), color = Colors.BLUE);
  # Experiment3.__plotRankFrequencies(xData = c(1, 2, 3, 4, 5), yData = c(147, 24, 7, 5, 5), color = Colors.VIOLET);
  
  # Plotter.createRanksHistograms(c(97, 21, 9, 7, 7), c(142, 24, 9, 5, 3),
  #                               c(141, 11, 11, 4, 1), c(147, 24, 7, 5, 5), titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
  # 
  # # length 15
  # Plotter.createRanksHistograms(c(125, 19, 8, 4, 6), c(161, 19, 4, 7, 4),
  #                               c(166, 9, 4, 3, 4), c(171, 8, 8, 5, 1),
  #                               titles = MultiTitles.RANK_HISTOGRAMS_3, Palettes.NATURE);
}


#########################################################
#' Computes histograms of the top years for each sample.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param pass {numerical} the pass which should be loaded
#' @param scoreType {string} the score-type which should be used for computation
#' @param topYearsCount {numerical} how many years considered for each profile
Experiment3.__computeYearFrequencyHistograms <- function(passPathGenericName, scoresFolderGenericName, 
                                                         pass, scoreType, topYearsCount) {
  passesPath <- paste(passPathGenericName, pass, sep = Symbols.EMPTY);
  
  scoresPerSample <- Loader.readInIndividualScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                         scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
  
  topYearsPerSample <- Analyzer.getTopYearsPerSample(scoresPerSample$scoresTables, topYearsCount);
  
  # iterate over each sample
  for (i in 1:length(topYearsPerSample)) {
    # retrieve
    fileName <- scoresPerSample$fileNames[[i]];
    plotFileName <- stringr::str_sub(fileName, start = 1, 
                                     end = nchar(fileName)-nchar(Strings.TO_REPLACE_1));
    plotPartName <- stringr::str_replace_all(plotFileName, Strings.TO_REPLACE_2, Strings.REPLACEMENT);
    plotFileName <- stringr::str_replace_all(plotFileName, Strings.TO_REPLACE_2, Strings.REPLACEMENT_2);
    
    sampleYears <- topYearsPerSample[[i]];
    
    # create
    title <- paste(Titles.YEAR_COUNT, plotPartName);
    Plotter.plotGeneralHistogram(unlist(sampleYears), title, 
                                 Strings.YEARS, Strings.FREQUENCY, TRUE, FALSE, plotFileName);
  }
}


#########################################################
#' Computes the box plots for the histogram based voting approach
#' out of the ranks from several computation passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param palette {vector} the colors palette which should be used
#' @param topYearsCount {numerical} how many years considered for each profile 
#' (warning: topYearsCount <= sampleLength has to hold in powerset approach)
#' @param powersetApproach {logical} tells if the powerset approach have to be used 
#' i.e. computing a score for every subset (hint: top-years count does not used anymore, if activated)
#' @param doubleWeighting {logical} tells if the years with ranks > 1 should be weighted down
#' @param logarithmicWeighting {logical} tells if the weights should be logarimized and rounded
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
Experiment3.__computeBoxPlots <- function(passPathGenericName, scoresFolderGenericName,
                                          passes, scoreTypes, palette, topYearsCount, 
                                          powersetApproach, doubleWeighting = FALSE,
                                          logarithmicWeighting = FALSE, minimumLength = -1) {
  
  combinedRanksPerScoreType <- Encoder.computeBoxDataPlot(passPathGenericName, scoresFolderGenericName,
                                                          passes, scoreTypes, palette, topYearsCount, 
                                                          powersetApproach, doubleWeighting,
                                                          logarithmicWeighting, minimumLength);
  
  # (D): plot scores
  data <- Plotter.plotBoxPlots(combinedRanksPerScoreType, scoreTypes, 
                               Titles.VIOLINE_RANKS_PLOT_3,
                               TRUE, FALSE, Extensions.PDF, palette);
}


#########################################################
#' Computes the histogram goodness by measuring 
#' the distance from best to second best.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} contains the scoreType
#' @param topYearsCount {numerical} how many years considered for each profile
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
#' @param xLimits {vector} the 2D vector specifying upper and lower bound
#' @param yLimits {vector} the 2D vector specifying upper and lower bound
#' @param print {logical} show the plot or not (default)
#' 
#' @return {list(xData, yData)} the x- and y-data from second plot
Experiment3.__evaluateHistogramGoodness <- function(passPathGenericName, scoresFolderGenericName,
                                                    passes, scoreType, topYearsCount, minimumLength,
                                                    xLimits, yLimits, print = TRUE) {
  
  data <- Encoder.createGoodnessHistogramEncoding(passPathGenericName, scoresFolderGenericName,
                                                  passes, scoreType, topYearsCount, 
                                                  minimumLength);
  
  # first evaluation
  xData <- data$xData;
  yData <- data$yData;
  yOutcast <- data$yOutcast;
  
    # encode in graph
  plot1 <- Plotter.plotGeneralScatterPlot(xData = xData, 
                                          yData = yData, 
                                          group = rep(Strings.POINT, length(yData)), 
                                          plotTitle = Titles.QUALITY_FACTOR, 
                                          xAxisLabel = Strings.DELTA_PEAK, 
                                          yAxisLabel = Strings.IS_RANK_1, 
                                          groupLabel = Strings.SAMPLES, 
                                          colors = Colors.VIOLET, 
                                          shapes = c(18, 15, 17, 18),
                                          print = FALSE, 
                                          save = FALSE);
  
  print(Strings.OUTLIERS);
  print(paste(Strings.MEAN_VALUE, mean(yOutcast)));
  print(paste(Strings.NUMBER, length(yOutcast)));
  print(paste(Strings.NUMBER_RANK_1, length(which(yOutcast == 1))));
  
  # second evaluation
    # create ranges
  value0to25 <- which(xData <= 25);
  value25to50 <- which(xData > 25 & xData <= 50);
  value50to75 <- which(xData > 50 & xData <= 75);
  value75to100 <- which(xData > 75 & xData <= 100);
  
  value100to125 <- which(xData > 100 & xData <= 125);
  value125to150 <- which(xData > 125 & xData <= 150);
  value150to175 <- which(xData > 150 & xData <= 175);
  value175to200 <- which(xData > 175 & xData <= 200);
  
  value200to225 <- which(xData > 200 & xData <= 225);
  
    # compute averages at these peaks
  average0to25 <- mean(yData[value0to25]);
  average25to50 <- mean(yData[value25to50]);
  average50to75 <- mean(yData[value50to75]);
  average75to100 <- mean(yData[value75to100]);
  
  average100to125 <- mean(yData[value100to125]);
  average125to150 <- mean(yData[value125to150]);
  average150to175 <- mean(yData[value150to175]);
  average175to200 <- mean(yData[value175to200]);
  
  average200to225 <- mean(yData[value200to225]);
  
    # encode
  xData <- c(12.5, 37.5, 62.5, 87.5, 
             112.5, 137.5, 162.5, 187.5, 
             212.5);
  yData <- c(average0to25, average25to50,
             average50to75, average75to100,
             average100to125, average125to150,
             average150to175, average175to200,
             average200to225);
  
  plot2 <- Plotter.plotGeneralLinePlot(xData = xData, 
                                       yData = yData, 
                                       plotTitle = Titles.MEANS_IN_RANGES, 
                                       xAxisLabel = Strings.DELTA_PEAK, 
                                       yAxisLabel = Strings.MEANS, 
                                       print = FALSE, 
                                       save = FALSE,
                                       points = TRUE,
                                       xLimits = xLimits,
                                       yLimits = yLimits,
                                       fitCurve = FALSE);
  
  if (print) {
    combinedPlot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, nrow = 2);
    print(combinedPlot);
  }
  
  return(list(xData = xData, yData = yData));
}


#########################################################
#' Plots the frequencies for the ranks.
#'
#' @param xData {vector} the data for the x-axis
#' @param yData {vector} the corresponding data for the y-axis
#' @param color {string} the color used for the curve
#' @export
Experiment3.__plotRankFrequencies <- function(xData, yData, color) {
  Plotter.plotGeneralLinePlot(xData, yData, plotTitle = Titles.WIDNOW_QUALITY, 
                              xAxisLabel = Strings.PEAK_RANK, yAxisLabel = Strings.FREQUENCY, TRUE, FALSE,
                              points = TRUE, xLimits = c(0.5, 5), yLimits = c(0, 170), percentage = FALSE, color = color);
}