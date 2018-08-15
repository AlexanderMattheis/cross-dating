#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_5_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes fith reproducible experiment which was used to create the Per-Tree Approach results.

#########################################################
#' Starts the per-tree approach.
Experiment5.start <- function() {
  # per tree approach
  # Experiment5.__storeComputedCounts(1, 1916, 2004);
  
  # Experiment5.__createHistogramPerSampleTrees(Paths.PASSES_BUCKET_MIN, Paths.PASSES_COUNTS, 
  #                                             Files.SCORES_FOLDER_NAME, list(1, 2, 3, 4, 5),  Defaults.TASK_D);
  
  # check runtime
  # runtimeDifferences <- list();
  # for (i in 1:5) {
  #   difference <- Experiment5.__computePerTreeScores(Defaults.TASK_C, 1);
  #   runtimeDifferences <- rlist::list.append(runtimeDifferences, difference);
  # }
  # 
  # for (i in seq_along(runtimeDifferences)) {
  #   print(runtimeDifferences[[i]]);
  # }
  
  # compute ranks
  Experiment1.createSummarizedViolinePlot(Paths.PASSES_PER_TREE, Files.SCORES_FOLDER_NAME,
                                          list(1),  list("d"),
                                          Palettes.NATURE);
  
  # check against alternative implementation
  # a <- Experiment5.computePerTreeScoresAlternative1(1);
  # b <- Experiment5.computePerTreeScoresAlternative2(1);
  # 
  # print("Test");
  
  # same approach with ring-widths
  # Experiment5.__computePerTreeRingWidthScores(Defaults.CORRELATION_T_VALUE, Defaults.TASK_T_VALUE_CHAR, 7);
  # Experiment5.__computePerTreeRingWidthScores(Defaults.CORRELATION_KENDALL, Defaults.TASK_KENDALL_CHAR, 7);
  # Experiment5.__computePerTreeRingWidthScores(Defaults.CORRELATION_PEARSON, Defaults.TASK_PEARSON_CHAR, 7);
  # Experiment5.__computePerTreeRingWidthScores(Defaults.CORRELATION_SPEARMAN, Defaults.TASK_SPEARMAN_CHAR, 7);
  
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_PER_TREE_MAX_DENSITY, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5),  list("p", "r", "t", "v"),
  #                                         Palettes.OTHER_NATURE, decreasing = TRUE);
  
  # length dependency best method
  # Experiment1.createViolineForMethod(list(Paths.PASSES_BUCKET_MIN_LENGTH_1, Paths.PASSES_PER_TREE_LENGTH_5,
  #                                         Paths.PASSES_PER_TREE, Paths.PASSES_PER_TREE_LENGTH_15), list("01", "05", "10", "15"),
  #                                    Files.SCORES_FOLDER_NAME, list(1, 2, 3, 4, 5, 6, 7), Defaults.TASK_C, Palettes.BLUE);
  
  # tested approach
  # Experiment5.__computePerTreeScoresTrivial(1);
  
  # Experiment5.__computePerTreeScoresWithSubdivision(pass = 7, scoreType = Defaults.TASK_A, numDivisions = 15);
  # Experiment5.__computePerTreeScoresWithSubdivision(pass = 7, scoreType = Defaults.TASK_B, numDivisions = 15);
  # Experiment5.__computePerTreeScoresWithSubdivision(pass = 7, scoreType = Defaults.TASK_C, numDivisions = 15);
  # Experiment5.__computePerTreeScoresWithSubdivision(pass = 7, scoreType = Defaults.TASK_D, numDivisions = 15);
  
  # Experiment5.__plotSplitDependency(c(1, 2, 5, 10), c(103, 114, 103, 98), Strings.QUALITY_2, Colors.GREEN);
  # Experiment5.__plotSplitDependency(c(1, 2, 5, 10), c(93, 126, 137, 144), Strings.QUALITY_2, Colors.RED);
  # Experiment5.__plotSplitDependency(c(1, 2, 5, 10), c(136, 147, 145, 140), Strings.QUALITY_2, Colors.BLUE);
  # Experiment5.__plotSplitDependency(c(1, 2, 5, 10), c(91, 119, 138, 147), Strings.QUALITY_2, Colors.VIOLET);
  
  # Experiment5.__plotSplitDependency(c(1, 3, 5, 15), c(133, 146, 156, 126), Strings.QUALITY_2, Colors.GREEN);
  # Experiment5.__plotSplitDependency(c(1, 3, 5, 15), c(126, 169, 179, 162), Strings.QUALITY_2, Colors.RED);
  # Experiment5.__plotSplitDependency(c(1, 3, 5, 15), c(159, 179, 183, 165), Strings.QUALITY_2, Colors.BLUE);
  # Experiment5.__plotSplitDependency(c(1, 3, 5, 15), c(143, 164, 174, 172), Strings.QUALITY_2, Colors.VIOLET);
}


#########################################################
#' Computes for all samples the counts for each position in method d) in the bucket-chronology
#' and stores this counts
#'
#' @param pass {numeric} the pass which should be loaded
#' @param startYear {numeric} the start year from which on the bucket-chronology should be analyzed
#' @param endYear {numeric} the last analyzed year in the bucket-chronology
Experiment5.__storeComputedCounts <- function(pass, startYear, endYear) {
  # retrieve path
  passesPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  perYearProfilesCleaned <- collectedData$profilesPerYear;
  perYearNamesCleaned <- collectedData$namesPerYear;
  
  countsPerSample <- Analyzer.computeCountsPerSample(perYearProfilesCleaned, perYearNamesCleaned, sampleData);
  
  Storer.storeCounts(countsPerSample, sampleData, 
                     startYear:endYear, 
                     Defaults.TASK_D, Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Analyzes from how many trees, the distance 
#' for the rank 1 rated samples is constructed.
#' Therefore, it looks on the number of trees
#' per year from which the minimum distance is constructed.
#'
#' @param passPathGenericName {string} the generic name of the scores pass path (without the pass number)
#' @param passPathGenericNameCounts {string} the generic name of the counts pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType for which the histogram should be computed
#' @export
Experiment5.__createHistogramPerSampleTrees <- function(passPathGenericName, passPathGenericNameCounts, 
                                                        scoresFolderGenericName, passes, scoreType) {
  # load scores
  counts <- Encoder.createHistogramPerSampleTreesEncoding(passPathGenericName, passPathGenericNameCounts, 
                                                          scoresFolderGenericName, passes, scoreType, 1);
  
  Plotter.plotGeneralHistogram(unlist(counts), Titles.TREES_PER_RANK_1, 
                               Strings.NUMBER_OF_TREES, Strings.FREQUENCY, 
                               TRUE, FALSE, Defaults.HISTOGRAM, FALSE, 
                               xLimits = c(0, 11), yLimits = c(0, 60));
}


#########################################################
#' Computes the per tree approach.
#' Instead of comparing profiles against buckets,
#' it is tested only against profiles from the same tree
#' to reduce the runtime.
#'
#' @param scoreType {list} the scoreType for which the histogram should be computed
#' @param pass {numerical} the pass which should be loaded
Experiment5.__computePerTreeScores <- function(scoreType, pass) {
  # retrieve path
  passesPath <- paste(Paths.PASSES_PER_TREE, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  
  # runtime optimizations #
  # There will be multiple separate chronologies, each with different start- & end-year.
  # And a test-sample should not be tested against the chronology if it is outside the chronology.
  # That is why the start- and end-years are stored per tree.
  yearsData <- Analyzer.getStartAndEndYears(cleanedCurvesData);
  startYears <- yearsData$startYears;
  endYears <- yearsData$endYears;
  
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  minStartYear <- collectedData$minStartYear;
  maxEndYear <- collectedData$maxStartYear;
  
  # search samples in bucket chronology
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, sampleData,
                                                        FALSE, FALSE, startYears, endYears, minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, sampleData,
                                                        FALSE, TRUE, startYears, endYears, minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, sampleData,
                                                        TRUE, FALSE, startYears, endYears, minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, sampleData,
                                                        TRUE, TRUE, startYears, endYears, minStartYear:maxEndYear);
  }
  
  print(Strings.FINISHED);
  
  Storer.storeIndividualScores(scoresPerSample, sampleData,
                               minStartYear:maxEndYear,
                               scoreType, Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Computes the per tree approach [for testing purposes].
#' Instead of comparing profiles against buckets,
#' it is tested only against profiles from the same tree
#' to reduce the runtime.
#'
#' @param pass {numerical} the pass which should be loaded
#' @param pathTrees {string} the path to the trees
#' @param pathDataSet {string} the path to the samples
#' @return {vector} the scores of the first sample
#' @export
Experiment5.computePerTreeScoresAlternative1 <- function(pass, pathTrees = Paths.PASSES_PER_TREE, 
                                                         pathDataSet = Paths.OSTALB_DATASET) {
  # retrieve path
  passesPath <- paste(pathTrees, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(pathDataSet, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  
  # runtime optimizations #
  # There will be multiple separate chronologies, each with different start- & end-year.
  # And a test-sample should not be tested against the chronology if it is outside the chronology.
  # That is why the start- and end-years are stored per tree.
  yearsData <- Analyzer.getStartAndEndYears(cleanedCurvesData);
  startYears <- yearsData$startYears;
  endYears <- yearsData$endYears;
  
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  minStartYear <- collectedData$minStartYear;
  maxEndYear <- collectedData$maxStartYear;
  
  sampleData <- list(sampleData[[1]]);
  
  # old <- Sys.time();
  scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, sampleData,
                                                      FALSE, FALSE, startYears, endYears, minStartYear:maxEndYear);
  # difference <- Sys.time() - old;
  # print(difference);
  
  print(Strings.FINISHED);
  
  return(unlist(scoresPerSample));
}


#########################################################
#' Computes the per tree approach [alternative implementation for testing purposes].
#' Instead of comparing profiles against buckets,
#' it is tested only against profiles from the same tree
#' to reduce the runtime.
#'
#' @param pass {numerical} the pass which should be loaded
#' @param pathTrees {string} the path to the trees
#' @param pathDataSet {string} the path to the samples
#' @return {vector} the scores of the first sample
#' @export
Experiment5.computePerTreeScoresAlternative2 <- function(pass, pathTrees = Paths.PASSES_PER_TREE,
                                                         pathDataSet = Paths.OSTALB_DATASET) {
  # retrieve path
  passesPath <- paste(pathTrees, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(pathDataSet, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology

  perYearScores <- new.env();

  sample <- samples[[1]][[1]];
  
  # old <- Sys.time();
  for (i in seq_along(cleanedCurvesData)) {
    chronology <- cleanedCurvesData[[i]];
    # print(chronology$name);
    perYearScores <- Analyzer.computeScoresOfSample(chronology, sample, perYearScores);  # equals method a)
  }
  
  # difference <- Sys.time() - old;
  # print(difference);
  
  # iterate over all per year scores to get final scores
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  minStartYear <- collectedData$minStartYear;
  maxEndYear <- collectedData$maxStartYear;
  
  minimaScores <- list();
  
  for (year in minStartYear:maxEndYear) {
    minimaScore <- Inf;
    
    if (!is.null(perYearScores[[toString(year)]])) {
      minimaScore <- min(unlist(perYearScores[[toString(year)]]));
    } 

    minimaScores <- rlist::list.append(minimaScores, minimaScore);
  }
  
  return(unlist(minimaScores));
}


#########################################################
#' Computes the per tree approach on ring-widths.
#'
#' @param method {string} the coefficient type with which the scores should be computed
#' @param scoreName {string} the name for the score
#' @param pass {numerical} the pass which should be loaded
Experiment5.__computePerTreeRingWidthScores <- function(method, scoreName, pass) {
  passesPath <- paste(Paths.PASSES_PER_TREE_MAX_DENSITY_LENGTH_15, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_MAXIMUM_DENSITIES, Symbols.REAL_DATA_SEPARATOR, TRUE, Extensions.META);
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_MAXIMUM_DENSITIES, Symbols.REAL_DATA_SEPARATOR, TRUE, 
                                    columnNameWidth = Defaults.DENSITY_NAME);
  samples <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                 Symbols.REAL_DATA_SEPARATOR);
    
  sampleData <- CurvesMiner.reconvertSamples(samples);  # convert in correct format
  sampleData <- CurvesMiner.extractCorrespondingWidths(curvesData, sampleData);
  
  cleanedCurvesData <- CurvesMiner.subtractSampleData(curvesData, sampleData, TRUE);  # train set
  
  # runtime optimizations #
  # There will be multiple separate chronologies, each with different start- & end-year.
  # And a test-sample should not be tested against the chronology if it is outside the chronology.
  # That is why the start- and end-years are stored per tree.
  yearsData <- Analyzer.getStartAndEndYears(cleanedCurvesData);
  startYears <- yearsData$startYears;
  endYears <- yearsData$endYears;
  
  collectedData <- CurvesMiner.collectWidthsData(cleanedCurvesData);
  minStartYear <- collectedData$minStartYear;
  maxEndYear <- collectedData$maxEndYear;
    
  # old <- Sys.time();
  scoresPerSample <- Analyzer.computeCoefficientsPerSample2(cleanedCurvesData, sampleData$samples, method, 
                                                            startYears, endYears, minStartYear:maxEndYear);
  # difference <- Sys.time() - old;
  # print(difference);
  
  Storer.storeScores(scoresPerSample, samples, 
                     minStartYear:maxEndYear, 
                     scoreName, Symbols.REAL_DATA_SEPARATOR); 
}


#########################################################
#' Computes the tree scores trivially
#' by going through each chronology separatly.
#' So with each chronology a score is generated for a year
#' and then the minimum score in each year is selected
#' as the final score. Just for testing purposes.
#' 
#' @param pathTrees {string} the path to the samples
#' @param pathDataSet {string} the path to the dataset with the trees
#' @param pass {numerical} the pass which should be loaded
Experiment5.__computePerTreeScoresTrivial <- function(pass, pathTrees = Paths.PASSES_PER_TREE,
                                                      pathDataSet = Paths.OSTALB_DATASET) {
  # retrieve path
  passesPath <- paste(pathTrees, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(pathDataSet, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create bucket chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  minStartYear <- collectedData$minStartYear;
  maxEndYear <- collectedData$maxStartYear;
  
  perSampleMinimaScores <- list();
  
  # iterate over all samples
  for (j in seq_along(samples[[1]])) {
    sample <- samples[[1]][[j]];
    
    perYearScores <- new.env();
    
    # iterate over all chronologies
    for (i in seq_along(cleanedCurvesData)) {
      chronology <- cleanedCurvesData[[i]];
      perYearScores <- Analyzer.computeScoresOfSample(chronology, sample, perYearScores);  # equals method a)
    }
    
    minimaScores <- list();
    
    for (year in minStartYear:maxEndYear) {
      minimaScore <- Inf;
      
      if (!is.null(perYearScores[[toString(year)]])) {
        minimaScore <- min(unlist(perYearScores[[toString(year)]]));
      } 
      
      minimaScores <- rlist::list.append(minimaScores, minimaScore);
    }
    
    perSampleMinimaScores <- rlist::list.append(perSampleMinimaScores, minimaScores);
  }
  
  # difference <- Sys.time() - old;
  # print(difference);
  
  View(perSampleMinimaScores);
}


#########################################################
#' Subdivides the samples into parts and then computes with the first part scores
#' and then with second part. The scores of both parts are then added together
#' such that you get a final score for a sample which was subdivided 
#' for example in to two parts that were separatly tested.
#' 
#' @param pass {numerical} the pass which should be loaded
#' @param scoreType {string} the score-type that should be used (a-d)
#' @param numDivisions {numeric} the number of divisions
Experiment5.__computePerTreeScoresWithSubdivision <- function(pass, scoreType, numDivisions) {
  # retrieve path
  passesPath <- paste(Paths.PASSES_PER_TREE_SUBDIV_15_LENGTH_15, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  # load data
  allCurves <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                    Symbols.REAL_DATA_SEPARATOR);
  
  # create chronology
  samples <- CurvesMiner.reconvertSamples(sampleData);  # convert in correct format
  cleanedCurvesData <- CurvesMiner.subtractSampleData(allCurves, samples);  # removing samples from data for chronology
  
  # runtime optimizations #
  # There will be multiple separate chronologies, each with different start- & end-year.
  # And a test-sample should not be tested against the chronology if it is outside the chronology.
  # That is why the start- and end-years are stored per tree.
  yearsData <- Analyzer.getStartAndEndYears(cleanedCurvesData);
  startYears <- yearsData$startYears;
  endYears <- yearsData$endYears;
  
  collectedData <- CurvesMiner.collectData(cleanedCurvesData);
  minStartYear <- as.numeric(collectedData$minStartYear);
  maxEndYear <- as.numeric(collectedData$maxStartYear);
  
  # subdivide sample data
  data <- CurvesMiner.subdivideSampleData(sampleData, numDivisions);
  
  subsampleLength <- data$subsampleLength;
  subdividedSampleData <- data$subdividedSampleData;
  
  divisions <- list();
  
  # iterate per subdivision -> hint: (minStartYear + (i-1)*subsampleLength) won't work due to optimizations!
  for (i in seq_along(subdividedSampleData)) {
    newSampleData <- subdividedSampleData[[i]];
    
    # search samples in bucket chronology
    if (scoreType == Defaults.TASK_A) {
      scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, newSampleData,
                                                          FALSE, FALSE, startYears, endYears, 
                                                          minStartYear:maxEndYear);
    } else if (scoreType == Defaults.TASK_B) {
      scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, newSampleData,
                                                          FALSE, TRUE, startYears, endYears,
                                                          minStartYear:maxEndYear);
    } else if (scoreType == Defaults.TASK_C) {
      scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, newSampleData,
                                                          TRUE, FALSE, startYears, endYears,
                                                          minStartYear:maxEndYear);
    } else if (scoreType == Defaults.TASK_D) {
      scoresPerSample <- Analyzer.computeScoresPerSample4(cleanedCurvesData, newSampleData,
                                                          TRUE, TRUE, startYears, endYears,
                                                          minStartYear:maxEndYear);
    }
    
    divisions <- rlist::list.append(divisions, scoresPerSample);
  }
  
  length <- length(minStartYear:maxEndYear);
  
  finalScoresPerSample <- list();
  
  # iterate over all samples
  for (j in 1:length(scoresPerSample)) {
    finalScoresOfSample <- 0;
    
    # iterate over all division and sum them up
    for (i in 1:length(divisions)) {
      scoresPerSample <- divisions[[i]];
      
      if (i > 1) {
        scoresOfSample <- unlist(scoresPerSample[[j]])[-seq(i*subsampleLength-subsampleLength)];
      } else {
        scoresOfSample <- unlist(scoresPerSample[[j]]);
      }
      
      newLength <- length(scoresOfSample);
      
      paddingEnd <- length-newLength;
      scoresOfSample <- c(scoresOfSample, rep(NA, paddingEnd));
      finalScoresOfSample <- finalScoresOfSample + scoresOfSample;
    }
    
    finalScoresPerSample <- rlist::list.append(finalScoresPerSample, as.list(finalScoresOfSample));
  }
  
  print(Strings.FINISHED);
  
  Storer.storeIndividualScores(finalScoresPerSample, sampleData,
                               minStartYear:maxEndYear,
                               scoreType, Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Creates a line plot showing the dependency 
#' between length and mean ranks for example.
#'
#' @param xData {vector} the data for the x-axis
#' @param yData {vector} the corresponding data for the y-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @export
Experiment5.__plotSplitDependency <- function(xData, yData, yAxisLabel, color) {
  Plotter.plotGeneralLinePlot(xData, yData, plotTitle = Titles.WIDNOW_QUALITY, 
                              xAxisLabel = Strings.NUMBER_OF_SUBSAMPLES, yAxisLabel = yAxisLabel, TRUE, FALSE,
                              points = TRUE, xLimits = c(0, 1, 3, 5, 15), yLimits = c(120, 185), percentage = FALSE, color = color,
                              xBreaks = 6);
}