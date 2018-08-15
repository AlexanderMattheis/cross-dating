#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_2_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes second reproducible experiment which was used to create the Bucket Approach results.

#########################################################
#' Starts the experiment.
Experiment2.start <- function() {
  # create scores
  # Experiment2.storeComputedScores(Defaults.TASK_A, 1, min, sum);  # the inner and outer function
  # Experiment2.storeComputedScores(Defaults.TASK_B, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_C, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);

  # compute and plot ranks with the help of the scores
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                                         Palettes.NATURE);
  
  # length 15
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN_LENGTH_1, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                                         Palettes.NATURE);
  
  # (E): create Venn-Diagrams
  # plot rank 1
  # Experiment1.createVennDiagrams(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME, 5, list(1));
  # Experiment1.createVennDiagrams(Paths.PASSES_BUCKET_MIN_LENGHT_1, Files.SCORES_FOLDER_NAME, 1, list(1));
  
  # plot ranks with worst indices
  # Experiment1.createVennDiagrams(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME, 4, list(-1:-10), TRUE);
  
  # (H): create back-to-back histogram
  # Experiment1.createBackToBackHistogram(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME, list(1:10), list(1, 2, 3, 4, 5));
  # Experiment1.createBackToBackHistogram(Paths.PASSES_BUCKET_MIN_LENGHT_1, Files.SCORES_FOLDER_NAME, list(1:90), list(1,2,3,4,5));
  
  # (J): create bar plot
  # Experiment2.createRatioBarPlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                1:10, list(1, 2, 3, 4, 5), Defaults.CURVE_OSTALB_START_YEAR,
  #                                Defaults.CURVE_OSTALB_END_YEAR);
  
  # (K): create delta-score/rank plot
  # Experiment2.__createDifferencePlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5),  Defaults.TASK_D,
  #                                    Colors.VIOLET, c(18, 15, 17, 18));

  # create weighted plot
  # Experiment2.__createDifferencePlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5),  Defaults.TASK_D,
  #                                    Colors.VIOLET, c(18, 15, 17, 18), FALSE);
  
  # Experiment2.__evaluateScoreGoodness(passPathGenericName = Paths.PASSES_BUCKET_MIN,
  #                                     scoresFolderGenericName = Files.SCORES_FOLDER_NAME,
  #                                     passes = list(1, 2, 3, 4, 5),
  #                                     scoreType = Defaults.TASK_D,
  #                                     years = as.list(1916:2004),
  #                                     xLimits = c(0, 8),
  #                                     yLimits = c(-1.1, 1.1),
  #                                     range = 8/9,
  #                                     incrementFrequency = 9);
  
  # Experiment2.__evaluateScoreGoodnessWithPValues(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                                list(1, 2, 3, 4, 5),  list("d"), 1,
  #                                                xLimits = c(0, 0.106),
  #                                                yLimits =  c(-1.1, 1.1),
  #                                                range = 11/900,
  #                                                incrementFrequency = 9,
  #                                                secondPlotTitle = Titles.MEANS_IN_RANGES_11_900,
  #                                                print = TRUE);
  
  # Experiment2.__evaluateScoreGoodnessWithPValues(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                                list(1, 2, 3, 4, 5),  list("d"), 1,
  #                                                xLimits = c(0, 0.081),
  #                                                yLimits =  c(-1.1, 1.1),
  #                                                print = TRUE,
  #                                                range = 0.009,
  #                                                incrementFrequency = 9,
  #                                                secondPlotTitle = Titles.MEANS_IN_RANGES_0.009,
  #                                                computeDeltaScores = TRUE);
  
  # Experiment2.__plotScores(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                          list(1),  list("d"));
  
  
  # Experiment2.__plotScores(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                          list(1),  list("d"), computeDeltaScores = TRUE);
  
  # plots the mean curves from buckets approach and this voting approach together
  # Experiment2.__plotMeanCurvesTogether();
  
  # clustering approach
  # Experiment2.__computeBestClustering(pass = 1,
  #                                     plot = TRUE);
  
  # Experiment2.storeScoresForClusteredChronology(Defaults.TASK_A, 1, min, sum);
  # Experiment2.storeScoresForClusteredChronology(Defaults.TASK_B, 5, min, sum);
  # Experiment2.storeScoresForClusteredChronology(Defaults.TASK_C, 5, min, sum);
  # Experiment2.storeScoresForClusteredChronology(Defaults.TASK_D, 5, min, sum);
}


#########################################################
#' Computes for all samples the scores for each position in the bucket-chronology
#' and stores this scores.
#'
#' @param scoreType {string} the score-type which should be used for computation
#' @param pass {numerical} the pass which should be loaded
#' @param func {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the position
#' @param func2 {function} the function which should be applied on the 
#' per sample computed scores to get a score for the position
Experiment2.storeComputedScores <- function(scoreType, pass, func, func2) {
  # retrieve path
  passesPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  perYearProfilesCleaned <- CurvesMiner.collectData(cleanedCurvesData)$profilesPerYear;
  
  # search samples in bucket chronology
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        FALSE, FALSE, func, func2);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        FALSE, TRUE, func, func2);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        TRUE, FALSE, func, func2);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        TRUE, TRUE, func, func2);
  }
  
  print(Strings.FINISHED);
  
  Storer.storeIndividualScores(scoresPerSample, sampleData, 
                               Defaults.CURVE_OSTALB_START_YEAR:Defaults.CURVE_OSTALB_END_YEAR, 
                               scoreType, Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Creates a bar plot for the ratio of best and worst ranked years.
#' 
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param rankIndicesToLookAt {vector} the indices to look at from positive and negative site
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param startYear {number} the year at which the histogram should start
#' @param endYear {number} the year at which the histogram should end
#' @export
Experiment2.createRatioBarPlot <- function(passPathGenericName, scoresFolderGenericName, 
                                           rankIndicesToLookAt, passes, startYear, endYear) {
  xData <- startYear:endYear;
  ratios <- Encoder.createRatioBarPlotData(passPathGenericName, scoresFolderGenericName, 
                                           rankIndicesToLookAt, passes, startYear, endYear);
  
  Plotter.plotGeneralBarPlot(xData, unlist(ratios), Titles.INCREMENTED_YEARS_RATIO, 
                             Strings.YEARS, Strings.RATIO, TRUE, FALSE);
}


#########################################################
#' Creates delta-score/rank plot. So an xy plot with the given axes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param palette {vector} the colors palette which should be used
#' @param shapes {vector} the shapes palette which should be used
#' @param weighted {logical} tells if the lowest score should be weighted with its delta score or not 
#'
#' @example 
#' x = delta score = PredictionRank2Score-PredictionRank1Score (because PredictionRank2Score > PredictionRank1Score), 
#' y = (real) rank
Experiment2.__createDifferencePlot <- function(passPathGenericName, scoresFolderGenericName, 
                                               passes, scoreType, palette, shapes, weighted = FALSE) {
  
  data <- Encoder.createDifferenceScatterPlotData(passPathGenericName, scoresFolderGenericName, 
                                                  passes, scoreType, palette, shapes, weighted);
  
  xData <- unlist(data$xData);
  yData <- unlist(data$yData);
  group <- unlist(data$group);
  
  orderOfScores <- order(xData, decreasing = TRUE);
  sortedRanks <- yData[orderOfScores];
  firstHundredRanks <- sortedRanks[1:100];
  
  print(Strings.NUMBER_RANK_1);
  print(length(firstHundredRanks[firstHundredRanks == 1]));
  
  Plotter.plotGeneralScatterPlot(xData, yData, group, 
                                 Titles.RANKS_FOR_DELTA_SCORES, 
                                 Defaults.DELTA_SCORE,
                                 Strings.PREDICTED_CORRECT_YEAR_RANKS, Strings.RANKS, palette, shapes,
                                 TRUE, FALSE, xLimits = c(0,8));
}


#########################################################
#' Evaluates the delta-scores as quality meausure.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param years {list} the years for which the diagram has to be created
#' @param xLimits {vector} the 2D vector specifying upper and lower bound
#' @param yLimits {vector} the 2D vector specifying upper and lower bound
#' @param print {logical} show the plot or not (default)
#'
#' @return {list(xData, yData)} the x- and y-data from second plot
Experiment2.__evaluateScoreGoodness <- function(passPathGenericName, 
                                                scoresFolderGenericName,
                                                passes, scoreType, years,
                                                xLimits, yLimits, range, incrementFrequency, 
                                                print = TRUE) {
  
  data <- Encoder.createGoodnessHistogramEncoding2(passPathGenericName, 
                                                   scoresFolderGenericName,
                                                   passes, scoreType, years);
  # first evaluation
  xData <- data$xData;
  yData <- data$yData;
  
    # encode in graph
  plot1 <- Plotter.plotGeneralScatterPlot(xData = xData, 
                                          yData = yData, 
                                          group = rep(Strings.POINT, length(yData)), 
                                          plotTitle = Titles.QUALITY_FACTOR_2, 
                                          xAxisLabel = Defaults.DELTA_SCORE, 
                                          yAxisLabel = Strings.IS_RANK_1, 
                                          groupLabel = Strings.SAMPLES, 
                                          colors = Colors.VIOLET, 
                                          shapes = c(18, 15, 17, 18),
                                          print = FALSE, 
                                          save = FALSE,
                                          xLimits = xLimits,
                                          yLimits = yLimits);
  
  # second evaluation
  data <- Encoder.getMeanRangeEncoding(xData, yData, range, incrementFrequency);
  
  xData <- data$xData;
  yData <- data$yData;
  
  plot2 <- Plotter.plotGeneralLinePlot(xData = xData, 
                                       yData = yData, 
                                       plotTitle = Titles.MEANS_IN_RANGES_0_89, 
                                       xAxisLabel = Defaults.DELTA_SCORE, 
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
#' Evaluates the p-values as quality meausure.
#' Computes per sample the p-values 
#' for the rank 1 predictions scores
#' to meausure the significance.
#' 
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param threshold {numeric} the threshold for which the p-value has to be computed
#' @param xLimits {vector} the 2D vector specifying upper and lower bound
#' @param yLimits {vector} the 2D vector specifying upper and lower bound
#' @param range {numeric} tells the range size
#' @param incrementFrequency {numeric} tells how often it should be incremented
#' @param secondPlotTitle {string} the title for the lower plot
#' @param computeDeltaScores {logical} tells if delta scores should be used to compute p-values
#' @param print {logical} tells if the plot should be printed
#'
#' @export
Experiment2.__evaluateScoreGoodnessWithPValues <- function(passPathGenericName, scoresFolderGenericName, 
                                                           passes, scoreTypes, threshold, xLimits, yLimits, 
                                                           range, incrementFrequency, secondPlotTitle, 
                                                           computeDeltaScores = FALSE, print) {
  
  data <- Encoder.createGoodnessHistogramEncoding3(passPathGenericName, scoresFolderGenericName, 
                                                   passes, scoreTypes, threshold, computeDeltaScores);
  
  # first evaluation
  xData <- data$xData;
  yData <- data$yData;
  
  # encode in graph
  plot1 <- Plotter.plotGeneralScatterPlot(xData = xData, 
                                          yData = yData, 
                                          group = rep(Strings.POINT, length(yData)), 
                                          plotTitle = Titles.QUALITY_FACTOR_3, 
                                          xAxisLabel = Defaults.P_VALUES, 
                                          yAxisLabel = Strings.IS_RANK_1, 
                                          groupLabel = Strings.SAMPLES, 
                                          colors = Colors.VIOLET, 
                                          shapes = c(18, 15, 17, 18),
                                          print = TRUE, 
                                          save = FALSE,
                                          xLimits = xLimits);
  
  # second evaluation
  data <- Encoder.getMeanRangeEncoding(xData, yData, range, incrementFrequency);
  
  xData <- data$xData;
  yData <- data$yData;
  
  plot2 <- Plotter.plotGeneralLinePlot(xData = xData, 
                                       yData = yData, 
                                       plotTitle = secondPlotTitle, 
                                       xAxisLabel = Defaults.P_VALUES, 
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
#' Plots score-distribution from several passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param computeDeltaScores {vector} tells if delta scores should be plotted instead of scores
#'
#' @export
Experiment2.__plotScores <- function(passPathGenericName, scoresFolderGenericName, 
                                     passes, scoreTypes, deltaScores, computeDeltaScores = FALSE) {
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);

    # iterate over each score type
    for (j in 1:length(scoreTypes)) {
      # get score-type
      scoreType <- scoreTypes[[j]];
      
      # read in score-type scores per sample
      scoresPerSampleData <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                       scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
      scoresPerSample <- scoresPerSampleData$scoresData;
      fileNamesPerSample <- scoresPerSampleData$fileNames;
      
      # iterate over all per sample scores
      for (k in 1:length(scoresPerSample)) {
        scoresOfSample <- unlist(scoresPerSample[[k]]);
        fileNameOfSample <- fileNamesPerSample[[k]];
        
        if (computeDeltaScores) {
          scoresOfSample <- Analyzer.getDeltaScores(scoresOfSample);
        }
          
        Plotter.plotGeneralHistogram(data = scoresOfSample,
                                     plotTitle = fileNameOfSample,
                                     xAxisLabel = Strings.SCORES,
                                     yAxisLabel = Strings.FREQUENCY,
                                     print = FALSE,
                                     save = TRUE,
                                     fileName = paste(fileNameOfSample, scoreType,
                                                      sep = Symbols.EMPTY),
                                     barLabels = FALSE,
                                     binWidth = 1,
                                     drawDensity = TRUE,
                                     percentage = FALSE);
      }
    }
  }
}


#########################################################
#' Creates goodness curves into the same plot for scores and p-values.
Experiment2.__plotMeanCurvesTogether <- function() {
  data1 <- Experiment2.__evaluateScoreGoodness(passPathGenericName = Paths.PASSES_BUCKET_MIN, 
                                               scoresFolderGenericName = Files.SCORES_FOLDER_NAME,
                                               passes = list(1, 2, 3, 4, 5), 
                                               scoreType = Defaults.TASK_D,
                                               years = as.list(1916:2004),
                                               xLimits = c(0, 8), 
                                               yLimits =  c(-1.1, 1.1),
                                               print = FALSE);
  
  xData1 <- (data1$xData - min(data1$xData)) / (max(data1$xData) - min(data1$xData))
  yData1 <- data1$yData;
  group1 <- rep(Strings.SCORES_MEAN, length(yData1));
  
  data2 <- Experiment2.__evaluateScoreGoodnessWithPValues(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
                                                          list(1, 2, 3, 4, 5),  list("d"), 1,
                                                          xLimits = c(0, 0.106),
                                                          yLimits =  c(-1.1, 1.1),
                                                          print = TRUE);
  
  xData2 <- (data2$xData - min(data2$xData)) / (max(data2$xData) - min(data2$xData));
  yData2 <- data2$yData;
  group2 <- rep(Strings.P_VALUES, length(yData2));
  
  # combine together
  xData <- c(xData1, xData2);
  yData <- c(yData1, yData2);
  group <- c(group1, group2);
  
  Plotter.plotMultipleLinePlots(xData = xData, 
                                yData = yData, 
                                group = group, 
                                plotTitle = Titles.TRUTH_MEANS, 
                                xAxisLabel = Strings.RELATIVE_QUALITY, 
                                yAxisLabel = Strings.MEANS, 
                                groupLabel = Strings.TRUTH,
                                print = TRUE, 
                                save = FALSE, 
                                labels = c(Strings.SCORES_MEAN, Strings.P_VALUES), 
                                colors = c(Palettes.VIOLET[[1]], Palettes.VIOLET[[4]]), 
                                yLimits <- c(-1.1, 1.1),
                                points = TRUE,
                                noLegend = TRUE);
}


#########################################################
#' Computes a clustering for the given data.
#'
#' @param pass {numerical} the pass which should be loaded
#' @param plot {logical} tells if a histogram should be created before and after clustering
Experiment2.__computeBestClustering <- function(pass, plot = FALSE) {
  # retrieve path
  passesPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  perYearProfilesCleanedData <- CurvesMiner.collectData(cleanedCurvesData);
  perYearProfilesCleaned <- perYearProfilesCleanedData$profilesPerYear;
  
  startYear <- perYearProfilesCleanedData$minStartYear;
  endYear <- perYearProfilesCleanedData$maxStartYear;
  
  # cluster cleaned per year profiles and then store them (to avoid a recomputation for each score type)
  afterClusteringProfiles <- Clustering.getPerYearProfilesPerCluster(perYearProfilesCleaned, startYear, 
                                                                     endYear, limit = 25, selectionNumber = 24);
  
  # store profiles from clusters
  Storer.storeProfilesEnvironment(afterClusteringProfiles, startYear:endYear);
  
  if (plot) {
    data <- Encoder.createHistogramsEncoding(startYear, endYear, perYearProfilesCleaned, afterClusteringProfiles);
    
    plot1 <- Plotter.plotGeneralHistogram(unlist(data$yearsBefore), Titles.PROFILES_PER_YEAR_BEFORE, 
                                          Strings.YEARS, Strings.FREQUENCY, TRUE, FALSE, 
                                          Defaults.HISTOGRAM, FALSE, yLimits = c(0, 60));
    
    plot2 <- Plotter.plotGeneralHistogram(unlist(data$yearsAfter), Titles.PROFILES_PER_YEAR_AFTER, 
                                          Strings.YEARS, Strings.FREQUENCY, TRUE, FALSE, 
                                          Defaults.HISTOGRAM, FALSE, yLimits = c(0, 60));
    
    combinedPlot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1, nrow = 2);
    print(combinedPlot);
  }
}


#########################################################
#' Computes for all samples the scores for each position in the 
#' chronology of buckets filled with profiles from clusters
#' and stores this scores.
#'
#' @param scoreType {string} the score-type which should be used for computation
#' @param pass {numerical} the pass which should be loaded
#' @param func {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the position
#' @param func2 {function} the function which should be applied on the 
#' per sample computed scores to get a score for the position
#' @export
Experiment2.storeScoresForClusteredChronology <- function(scoreType, pass, func, func2) {
  # retrieve path
  passesPath <- paste(Paths.PASSES_CLUSTERED, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # load clustered bucket chronology
  perYearProfilesData <- Loader.loadBuckets(paste(passesPath, Paths.BUCKETS, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
  
  perYearProfilesCleaned <- perYearProfilesData$perYearProfilesCleaned;
  startYear <- perYearProfilesData$startYear;
  endYear <- perYearProfilesData$endYear;
  
  # search samples in bucket chronology
  if (scoreType == Defaults.TASK_A) {
    old <- Sys.time();
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        FALSE, FALSE, func, func2);
    difference <- Sys.time() - old;
    print(difference);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        FALSE, TRUE, func, func2);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        TRUE, FALSE, func, func2);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfilesCleaned, sampleData,
                                                        TRUE, TRUE, func, func2);
  }
  
  Storer.storeIndividualScores(scoresPerSample, sampleData, 
                               startYear:endYear, 
                               scoreType, Symbols.REAL_DATA_SEPARATOR);
}