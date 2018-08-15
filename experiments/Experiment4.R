#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_4_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes fourth reproducible experiment which was used to create the Two-Step Approach results.

#########################################################
#' Starts the ring width based cross-dating against density-based cross-dating.
Experiment4.start <- function() {
  # stores samples and ring width consensus
  # Experiment4.__storeConsensusAndSamples(1);
  
  # compute correlation coefficients
  # Experiment4.__storeComputedScores(7, Defaults.TASK_T_VALUE_CHAR, Defaults.CORRELATION_T_VALUE);
  # Experiment4.__storeComputedScores(5, Defaults.TASK_GLEICHLAUF_CHAR, Defaults.CORRELATION_GLEICHLAUF);
  # Experiment4.__storeComputedScores(1, Defaults.TASK_SPEARMAN_CHAR, Defaults.CORRELATION_SPEARMAN);
  # Experiment4.__storeComputedScores(7, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  # Experiment4.__storeComputedScores(7, Defaults.TASK_KENDALL_CHAR, Defaults.CORRELATION_KENDALL);
  
  # compute bucket-based correlation coefficients
  # Experiment4.__storeComputedScores(5, Defaults.TASK_RING_WIDTH_BUCKET, Defaults.CORRELATION_SPEARMAN, TRUE);
  
  # compute ranks
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_MAX_DENSITY, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5),  list(Defaults.TASK_T_VALUE_CHAR,
  #                                                                    Defaults.TASK_PEARSON_CHAR,
  #                                                                    Defaults.TASK_SPEARMAN_CHAR,
  #                                                                    Defaults.TASK_KENDALL_CHAR),
  #                                         c(Palettes.OTHER_NATURE[4],
  #                                           Palettes.OTHER_NATURE[1],
  #                                           Palettes.OTHER_NATURE[2],
  #                                           Palettes.OTHER_NATURE[3]), TRUE, FALSE);
  
  # compute rank enhancement plot for the best e.g. 100 samples from the ring-width based approach
  # Experiment4.__computeRankEnhancementPlot(100, Paths.PASSES_MAX_DENSITY, Files.SCORES_FOLDER_NAME,
  #                                          list(1, 2, 3, 4, 5),  Defaults.TASK_PEARSON_CHAR);

  # extracts top years from ring-width based approach 
  # and looks up scores for these years in buckets-based approach
  # and create a rank-histogram for the correct years
  # Experiment4.__createsRankHistogramForTopYears(20, Paths.PASSES_MAX_DENSITY,
  #                                               Files.SCORES_FOLDER_NAME, list(1, 2, 3, 4, 5),
  #                                               Defaults.TASK_PEARSON_CHAR,
  #                                               Defaults.TASK_D);

  # pointer-year heuristic
  # Experiment4.__computePlotWithPointerYearsHeuristic(consensiPath = Paths.INPUT, 
  #                                                    consensusName = Files.OSTALB_CONSENSUS,
  #                                                    correlationCoefficientLimit = 0.96,
  #                                                    correlationCoefficientLimitForEvents = 0.83,
  #                                                    passPathGenericName = Paths.PASSES_BUCKET_MIN, 
  #                                                    scoresFolderGenericName = Files.SCORES_FOLDER_NAME, 
  #                                                    passes = list(1, 2, 3, 4, 5), 
  #                                                    scoreType = Defaults.TASK_D);
  
  # length dependant
  # Experiment1.createViolineForMethod(list(Paths.PASSES_MAX_DENSITY_LENGTH_5,
  #                                         Paths.PASSES_MAX_DENSITY, Paths.PASSES_MAX_DENSITY_LENGTH_15),
  #                                    list("05", "10", "15"),
  #                                    Files.SCORES_FOLDER_NAME, list(1, 2, 3, 4, 5, 6, 7),
  #                                    Defaults.TASK_PEARSON_CHAR, Palettes.RED);
  
  # compute max-density data
  # Experiment4.__storeMaxDensitiesOfProfiles();
  # Experiment4.__storeConsensusAndSamples(7);
}


#########################################################
#' Creates samples, a consensus and a consensus free of samples. 
#' Stores both on the disk.
#'
#' @param pass {numerical} the pass which should be loaded
Experiment4.__storeConsensusAndSamples <- function(pass) {
  # read in data
  passesPath <- paste(Paths.PASSES_LENGTH_15, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_WIDTHS, Symbols.REAL_DATA_SEPARATOR, TRUE, Extensions.META);
  # curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_MAXIMUM_DENSITIES, Symbols.REAL_DATA_SEPARATOR,
  #                                   columnNameWidth = Defaults.DENSITY_NAME, TRUE);
  
  samples <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                 Symbols.REAL_DATA_SEPARATOR);
  
  sampleData <- CurvesMiner.reconvertSamples(samples);  # convert in correct format
  sampleData <- CurvesMiner.extractCorrespondingWidths(curvesData, sampleData);
  
  # subtract sample data from curve data
  cleanedCurvesData <- CurvesMiner.subtractSampleData(curvesData, sampleData, TRUE);  # train set
  perYearWidthsCleaned <- CurvesMiner.collectWidthsData(cleanedCurvesData);
  masterChronologyCleaned <- CurvesMiner.computeCharacteristicsChronology(perYearWidthsCleaned$widthsPerYear,
                                                                          Defaults.CURVE_OSTALB_START_YEAR,
                                                                          Defaults.CURVE_OSTALB_END_YEAR);
  # store data
  Storer.storeConsecutiveSamples(sampleData, Symbols.REAL_DATA_SEPARATOR, TRUE);
  
  Storer.storeWidths(masterChronologyCleaned,
                     Defaults.MASTER_CHRONOLOGY_OSTALB_FILENAME,
                     Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Computes for all samples the scores for each position in the consensus
#' and stores this scores.
#' 
#' @param pass {numerical} the pass which should be loaded
#' @param scroreType {string} the score type name
#' @param method {string} the correlation method which should be used {"spearman", "pearson", "kendall"}
#' @param minDistanceApproach {logical} tells if the bucket based approach should be applied on the ring-width approach
Experiment4.__storeComputedScores <- function(pass, scoreType, method, minDistanceApproach = FALSE) {
  consensiPath <- paste(Paths.PASSES_MAX_DENSITY_LENGTH_15, pass, Paths.PATH, sep = Symbols.EMPTY);
  testSetPath <- paste(consensiPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY);
  passesPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  if (minDistanceApproach) {      # create bucket chronology
    curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_WIDTHS, Symbols.REAL_DATA_SEPARATOR, TRUE);
    samples <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY), 
                                   Symbols.REAL_DATA_SEPARATOR);

    sampleData <- CurvesMiner.reconvertSamples(samples);  # convert in correct format
    sampleData <- CurvesMiner.extractCorrespondingWidths(curvesData, sampleData);
    
    cleanedCurvesData <- CurvesMiner.subtractSampleData(curvesData, sampleData, TRUE);  # train set
    perYearWidthsCleaned <- CurvesMiner.collectWidthsData(cleanedCurvesData)$widthsPerYear;
    
    # old <- Sys.time();
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(perYearWidthsCleaned, sampleData$samples, method, minDistanceApproach);
    # difference <- Sys.time() - old;
    # print(difference);
    
    Storer.storeScores(scoresPerSample, samples, 
                       Defaults.CURVE_OSTALB_START_YEAR:Defaults.CURVE_OSTALB_END_YEAR, 
                       scoreType, Symbols.REAL_DATA_SEPARATOR); 
  } else {
    sampleData <- Loader.readInCurves(testSetPath, Symbols.REAL_DATA_SEPARATOR, TRUE, Extensions.CSV);
    masterChronologyData <- Loader.readInWidths(consensiPath, Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR, Extensions.CSV);

    # old <- Sys.time();
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(masterChronologyData$widths, sampleData, method, minDistanceApproach);
    # difference <- Sys.time() - old;
    # return(difference);
    # 
    Storer.storeScores(scoresPerSample, sampleData, 
                       Defaults.CURVE_OSTALB_START_YEAR:Defaults.CURVE_OSTALB_END_YEAR, 
                       scoreType, Symbols.REAL_DATA_SEPARATOR); 
  }
}


#########################################################
#' Computes a plot that shows the enhancement in the ranks
#' by using the buckets-based approach after the ring-width based
#' approach on a number of best ranked samples.
#'
#' @param number {numerical} the number of samples on which you want to look on
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the score type name of the ring-width based method
Experiment4.__computeRankEnhancementPlot <- function(number, passPathGenericName, 
                                                     scoresFolderGenericName, 
                                                     passes, scoreType) {
  
  data <- Encoder.createRankEnhancementPlotData(number, passPathGenericName, 
                                                scoresFolderGenericName, 
                                                passes, scoreType);
  
  x <- data$x;
  y <- data$y;
  g <- data$g;
  
  # create scatter plot
  Plotter.plotGeneralScatterPlot(x, unlist(y), unlist(g), Titles.ENHANCEMENT, 
                                 Strings.SAMPLES, Strings.PREDICTED_RANK, Titles.APPROACH, 
                                 c(Colors.VIOLET, Colors.PINK), c(18, 16),
                                 TRUE, FALSE, c(0.5, 1));
}


#########################################################
#' Creates a rank histogram after applying the buckets-based
#' approach on the top years from the ring-width based approach.
#'
#' @param topYears {numerical} the number of top scores that have to be retrieved 
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType of the ring-width approach
#' @param histogramScoreType {list} the scoreType of the buckets-based method
Experiment4.__createsRankHistogramForTopYears <- function(topYears, passPathGenericName, 
                                                          scoresFolderGenericName, passes, 
                                                          scoreType, histogramScoreType) {
  
  yData <- Encoder.createsRankHistogramForTopYears(topYears, passPathGenericName, 
                                                   scoresFolderGenericName, passes, 
                                                   scoreType, histogramScoreType);
  
  Plotter.plotGeneralHistogram(yData, Titles.RANKS, Strings.PREDICTED_CORRECT_YEAR_RANKS, Strings.FREQUENCY, 
                               TRUE, FALSE, Files.RANK_HISTOGRAM, TRUE, c(0, 200), xLimits = c(0, 20), color = Colors.VIOLET);
}


#########################################################
#' Computes a rank histogram under usage of the pointer-years heuristic.
#' First the years in the consensus-based approach determined 
#' which are potential pointer-year candidates.
#' Afterwards it is looked which samples contain event-years.
#' These samples are then searched only within the determined pointer-years.
#' All remaining samples are searched as before with the buckets approach.
#'
#' @param consensiPath {string} the path to the consensi
#' @param consensusName {string} the name of the consensus
#' @param correlationCoefficientLimit {numerical} the limit below which all years defined as pointer-years
#' @param correlationCoefficientLimitForEvents {numerical} the limit below which years defined as event-years
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType of the buckets-based method
Experiment4.__computePlotWithPointerYearsHeuristic <- function(consensiPath, consensusName,
                                                               correlationCoefficientLimit,
                                                               correlationCoefficientLimitForEvents,
                                                               passPathGenericName, 
                                                               scoresFolderGenericName, passes, 
                                                               scoreType) {
  
  encodingData <- Encoder.createRanksPointerHeuristicHistogramEncoding(consensiPath, consensusName,
                                                                       correlationCoefficientLimit,
                                                                       correlationCoefficientLimitForEvents,
                                                                       passPathGenericName, 
                                                                       scoresFolderGenericName, passes, 
                                                                       scoreType);
  yData <- encodingData$yData;
  nonCharacteristic <- encodingData$nonCharacteristicRanks;
  characteristic <- encodingData$characteristicRanks;
  outcasts <- encodingData$outcasts;
  
  Plotter.plotGeneralHistogram(yData, Titles.RANKS, Strings.RANKS, Strings.FREQUENCY, 
                               TRUE, FALSE, Files.RANK_HISTOGRAM, TRUE);
  
  print(Titles.CHARACTERISTIC);
  print(Strings.COUNT);
  print(length(nonCharacteristic));
  print(Strings.MEAN);
  print(mean(nonCharacteristic));
  
  print(Strings.SEPERATION);
  
  print(Strings.COUNT);
  print(length(characteristic));
  print(Strings.MEAN)
  print(mean(characteristic));
  
  print(Strings.OUTLIERS);
  print(outcasts);
}


#########################################################
#' Loads profiles and stores there maximum density.
Experiment4.__storeMaxDensitiesOfProfiles <- function() {
  # load profiles
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  
  maximumDensities <- CurvesMiner.getMaximumDensities(curvesData);
  
  # store the maxima
  Storer.storeMaxDensities(maximumDensities, Symbols.REAL_DATA_SEPARATOR);
}