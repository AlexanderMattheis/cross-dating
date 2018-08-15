#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_1_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes first reproducible experiment which was used to create the Consensus Approach results.

#########################################################
#' Starts the experiment.
Experiment1.start <- function() {
  # (A)/(F): extracts samples e.g. from Ostalb curves
  # (B): generates master-chronology with/without the samples
  # Experiment1.__storeConsensiAndSamples(FALSE);
  # Experiment1.__storeSubdividedSamples(Paths.PASSES, 5, 2);
  
  # (B): plot master-chronology with and without samples
  # Experiment1.__loadAndPlotConsensi(1);
  
  # (C): create rank violine scores
  # Experiment1.__storeComputedScores(Defaults.TASK_A, 1, 0);
  # Experiment1.__storeComputedScores(Defaults.TASK_B, 5, 0);
  # Experiment1.__storeComputedScores(Defaults.TASK_C, 1, 0);
  # Experiment1.__storeComputedScores(Defaults.TASK_D, 1, 0);
  # Experiment1.__storeComputedMicaScores(Defaults.TASK_E, 1);
  # Experiment1.__storeComputedMicaScores(Defaults.TASK_F, 1);
  # Experiment1.__storeComputedMicaScores(Defaults.TASK_G, 1);
  # Experiment1.__storeComputedMicaScores(Defaults.TASK_H, 1);

  # (D): compute and plot ranks with the help of the scores
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"),
  #                                         Palettes.NATURE);

  # Experiment1.createSummarizedViolinePlot(Paths.PASSES, Files.SCORES_FOLDER_NAME,
  #                                         list(1,2,3,4,5),  list("e", "f", "g", "h"),
  #                                         Palettes.NATURE);
  
  # (E): create Venn-Diagrams
  # plot rank 1
  # Experiment1.createVennDiagrams(Paths.PASSES, Files.SCORES_FOLDER_NAME, 1, list(1));
  
  # plot ranks with worst indices
  # Experiment1.createVennDiagrams(Paths.PASSES, Files.SCORES_FOLDER_NAME, 5, list(-1:-10), TRUE);
  
  # (F): create violine plots for different series-lengths
  # Experiment1.createViolineForMethod(list(Paths.PASSES_LENGTH_1, Paths.PASSES_LENGTH_5,
  #                                         Paths.PASSES, Paths.PASSES_LENGTH_15), list("01", "05", "10", "15"),
  #                                    Files.SCORES_FOLDER_NAME, list(1, 2, 3, 4, 5, 6, 7), Defaults.TASK_D, Palettes.VIOLET);
  
  # (G): create ROC-plots
  # Experiment1.createROC(Paths.PASSES, Files.SCORES_FOLDER_NAME,
  #                       list(1, 2, 3, 4, 5),  list("a", "b", "c", "d"));
  
  # (H): create back-to-back histogram
  # Experiment1.createBackToBackHistogram(Paths.PASSES, Files.SCORES_FOLDER_NAME, list(1:10), list(1, 2, 3, 4, 5));
  # Experiment1.createBackToBackHistogram(Paths.PASSES_LENGTH_1, Files.SCORES_FOLDER_NAME, list(1:90), list(1,2,3,4,5));
  
  # (I): create ratio plot
  # Experiment2.createRatioBarPlot(Paths.PASSES_LENGTH_1, Files.SCORES_FOLDER_NAME,
  #                                1:90, list(1, 2, 3, 4, 5), Defaults.CURVE_OSTALB_START_YEAR,
  #                                Defaults.CURVE_OSTALB_END_YEAR);
}


#########################################################
#' Creates samples, a consensus and a consensus free of samples. 
#' Stores both on the disk.
#'
#' @param computeReferenceChronology {logical} tells if a reference chronology has to be computed for comparisons
Experiment1.__storeConsensiAndSamples <- function(computeReferenceChronology) {
  # preprocessing
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);

  # (A): extract 50 samples or more from curvesData of lenght 10
  sampleData <- CurvesMiner.extractUniformlyConsecutiveSamples(numberOfSamples = 50, sampleLength = 10, curvesData);  # test set
  
  # (B): generate master-chronology without/with the samples from sampleData
  if (computeReferenceChronology) {
    perYearProfiles <- CurvesMiner.collectData(curvesData);
    masterChronology <- CurvesMiner.computeChronology(perYearProfiles$profilesPerYear,
                                                      Defaults.CURVE_OSTALB_START_YEAR,
                                                      Defaults.CURVE_OSTALB_END_YEAR, TRUE, c(0,2));
    
    # store data
    Storer.storeCurve2(masterChronology,
                       Defaults.CONSENSUS_CURVE_OSTALB_FILENAME,
                       Symbols.REAL_DATA_SEPARATOR);
    
    Storer.storeList(masterChronology$numberOfProfiles,
                     Files.OSTALB_CONSENSUS_NUM_PROFILES,
                     Symbols.REAL_DATA_SEPARATOR);
  }
  
  # subtract sample data from curve data
  cleanedCurvesData <- CurvesMiner.subtractSampleData(curvesData, sampleData);  # train set
  perYearProfilesCleaned <- CurvesMiner.collectData(cleanedCurvesData);
  masterChronologyCleaned <- CurvesMiner.computeChronology(perYearProfilesCleaned$profilesPerYear,
                                                           Defaults.CURVE_OSTALB_START_YEAR,
                                                           Defaults.CURVE_OSTALB_END_YEAR, FALSE);

  # store data
  Storer.storeConsecutiveSamples(sampleData, Symbols.REAL_DATA_SEPARATOR);

  Storer.storeCurve2(masterChronologyCleaned,
                     Defaults.MASTER_CHRONOLOGY_OSTALB_FILENAME,
                     Symbols.REAL_DATA_SEPARATOR);

  Storer.storeList(masterChronologyCleaned$numberOfProfiles,
                   Files.OSTALB_CONSENSUS_SAMPLE_FREE_NUM_PROFILES,
                   Symbols.REAL_DATA_SEPARATOR);

  Experiment1.__correctnessCheck1(curvesData, sampleData);
}


#########################################################
#' Loads samples from hard disk and stores them as subdivided samples again.
#' 
#' @param seriesLengthPath {string} the path to the series length
#' @param pass {numerical} the pass to look at
#' @param divisionFactor {numerical} the factor by which the sample length should be divided to create smaller samples
Experiment1.__storeSubdividedSamples <- function(seriesLengthPath, pass, divisionFactor) {
  passesPath <- paste(seriesLengthPath, pass, sep = Symbols.EMPTY);
  
  sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH),
                                    Symbols.REAL_DATA_SEPARATOR);
  
  if (divisionFactor > 0) sampleData <- CurvesMiner.subdivideSamples(sampleData, divisionFactor);
  
  Storer.storeConsecutiveSamples2(sampleData, Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Shows the user some information which can be used to see if
#' the data was preocessed correctly.
#'
#' @param curvesData {list} first parameter whose length should be displayed
#' @param sampleData {list(samples, years)} list of extracted samples
Experiment1.__correctnessCheck1 <- function(curvesData, sampleData) {
  numberOfCurves <- paste(Check.NUMBER_OF_CURVES, length(curvesData));
  numberOfSamples <- paste(Check.NUMBER_OF_SAMPLES, length(sampleData$samples));
  sampleLength <- paste(Check.SAMPLE_LENGTH, length(sampleData$samples[[1]]));
  
  print(numberOfCurves);
  print(numberOfSamples);
  print(sampleLength)
}


#########################################################
#' Loads the consensus, the sample-free consensus
#' and stores them year/profile-wise as plots on the hard disk.
#'
#' @param pass {numerical} the pass which should be loaded
Experiment1.__loadAndPlotConsensi <- function(pass) {
  # (B): plot master-chronology and consensus curve
  # load
  consensiPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  
  masterChronology <- Loader.readInProfiles(consensiPath, Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
  consensus <- Loader.readInProfiles(consensiPath, Files.OSTALB_CONSENSUS, Symbols.REAL_DATA_SEPARATOR);
  
  masterChronologyProfilesPerYear <- Loader.readInList(paste(consensiPath, Files.OSTALB_CONSENSUS_SAMPLE_FREE_NUM_PROFILES, 
                                                             Extensions.CSV, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
  consensusProfilesPerYear <- Loader.readInList(paste(consensiPath, Files.OSTALB_CONSENSUS_NUM_PROFILES, Extensions.CSV,
                                                      sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
  
  profilesNumPerYear <- list(masterChronology = masterChronologyProfilesPerYear, 
                             consensus = consensusProfilesPerYear);
  
  # plot
  Plotter.plotProfiles(consensus, masterChronology, Titles.WITH_AND_WITHOUT_SAMPLES_START,
                       Titles.CURVES, Strings.CONSENSUS_1_NAME, Strings.CONSENSUS_2_NAME,
                       TRUE, TRUE, Extensions.PDF, profilesNumPerYear, yLimits = c(0, 1.5));
}


#########################################################
#' Computes for all samples the scores for each position in the consensus
#' and stores this scores.
#'
#' @param scoreType {string} the score-type which should be used for computation
#' @param pass {numerical} the pass which should be loaded
#' @param divisionFactor {numerical} the factor by which the sample length should be divided to create smaller samples
Experiment1.__storeComputedScores <- function(scoreType, pass, divisionFactor) {
  print(scoreType);
  print(Titles.LOADING_FROM_DISK);
  
  consensiPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  testSetPath <- paste(consensiPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY);
  
  sampleData <- Loader.readInCurves(testSetPath, Symbols.REAL_DATA_SEPARATOR);
  masterChronologyData <- Loader.readInProfiles(consensiPath, Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);

  scoresPerSample <- list();
  
  if (divisionFactor > 0) sampleData <- CurvesMiner.subdivideSamples(sampleData, divisionFactor);
  
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, sampleData,
                                                       FALSE, FALSE);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, sampleData,
                                                       FALSE, TRUE);
  } else if (scoreType == Defaults.TASK_C) {
    # old <- Sys.time();
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, sampleData,
                                                       TRUE, FALSE);
    # difference <- Sys.time() - old;
    # return(difference);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, sampleData,
                                                       TRUE, TRUE);
  }
  
  Storer.storeScores(scoresPerSample, sampleData, 
                     masterChronologyData$years, 
                     scoreType, Symbols.REAL_DATA_SEPARATOR); 
}


#########################################################
#' Computes for all samples the MICA-scores for each position in the consensus
#' and stores this scores.
#' 
#' @param scoreType {string} the score-type which should be used for computation
#' @param pass {numerical} the pass which should be loaded
Experiment1.__storeComputedMicaScores <- function(scoreType, pass) {
  print(scoreType);
  print(Titles.LOADING_FROM_DISK);
  
  consensiPath <- paste(Paths.PASSES, pass, Paths.PATH, sep = Symbols.EMPTY);
  testSetPath <- paste(consensiPath, Paths.OSTALB_TESTSET_FOLDER, sep = Symbols.EMPTY);
  
  sampleData <- Loader.readInCurves(testSetPath, Symbols.REAL_DATA_SEPARATOR);
  masterChronologyData <- Loader.readInProfiles(consensiPath, Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
  
  scoresPerSample <- list();
  
  if (scoreType == Defaults.TASK_E) {
    scoresPerSample <- Analyzer.computeMicaScoresPerSample(masterChronologyData, sampleData,
                                                           FALSE, FALSE);
  } else if (scoreType == Defaults.TASK_F) {
    scoresPerSample <- Analyzer.computeMicaScoresPerSample(masterChronologyData, sampleData,
                                                           FALSE, TRUE);
  } else if (scoreType == Defaults.TASK_G) {
    scoresPerSample <- Analyzer.computeMicaScoresPerSample(masterChronologyData, sampleData,
                                                           TRUE, FALSE);
  } else if (scoreType == Defaults.TASK_H) {
    scoresPerSample <- Analyzer.computeMicaScoresPerSample(masterChronologyData, sampleData,
                                                           TRUE, TRUE);
  }
  
  Storer.storeScores(scoresPerSample, sampleData, 
                     masterChronologyData$years, 
                     scoreType, Symbols.REAL_DATA_SEPARATOR); 
}


#########################################################
#' Creates a violine plot out of the ranks from several computation passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param palette {vector} the colors palette which should be used
#' @param decreasing {logical} tells if the highest value should have lowest rank
#' @param exception {logical} tells if in the first round the decreasing parameter should be FALSE
#' @export
Experiment1.createSummarizedViolinePlot <- function(passPathGenericName, scoresFolderGenericName, 
                                                    passes, scoreTypes, palette, decreasing = FALSE,
                                                    exception = FALSE) {

  combinedRanksPerScoreType <- Encoder.createViolinePlotData(passPathGenericName, scoresFolderGenericName, 
                                                             passes, scoreTypes, decreasing,
                                                             exception); 
  # (D): plot ranks
  data <- Plotter.plotViolinePlots(combinedRanksPerScoreType, scoreTypes, Titles.VIOLINE_RANKS_PLOT, 
                                   TRUE, FALSE, Extensions.PDF, palette);
}


#########################################################
#' Creates Venn diagrams for the ranks of 4 measurement types.
#' So it is looked if the sets of samples for example for rank 1 have an overlap.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param pass {numerical} the pass which has to be used for computation
#' @param ranksToLookAt {list} the ranks for which the Venn Diagrams should be created 
#' (negative ranks indices gives you worst ranks i.e. -1 = last rank)
#' @param lookOnIndices {logical} tells if to look on indices or ranks i.e. changes ranksToLookAt
#' @export
Experiment1.createVennDiagrams <- function(passPathGenericName, scoresFolderGenericName, 
                                           pass, ranksToLookAt, lookOnIndices = FALSE) {
  
  data <- Encoder.createVennDiagramData(passPathGenericName, scoresFolderGenericName, 
                                        pass, ranksToLookAt, lookOnIndices);
  
  ranksPerScoreType <- data$ranksPerScoreType;
  fileNames <- data$fileNames;
  
  if (lookOnIndices) {  # look on indices of sorted ranks or ranks iterself
    Plotter.plotVennDiagram2(ranksPerScoreType, fileNames, unlist(ranksToLookAt), pass, TRUE);
  } else {
    Plotter.plotVennDiagram(ranksPerScoreType, fileNames, unlist(ranksToLookAt), pass, TRUE);
  }
}


#########################################################
#' Creates the violine plots for different series lengths.
#'
#' @param seriesLengthPaths {list} the paths to the different series lengths
#' @param lengths {list} the lengths of that paths
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' Hint: It can be be specified more than available.
#' @param scoreType {string} the score type at which to look at
#' @param palette {vector} the colors palette which should be used
#' @export
Experiment1.createViolineForMethod <- function(seriesLengthPaths, lengths, passPathGenericName, 
                                               passes, scoreType, palette) {

  ranksPerLength <- Encoder.createSeriesLengthViolinePlotData(seriesLengthPaths, lengths, passPathGenericName, 
                                                              passes, scoreType, palette);
  
  data <- Plotter.plotViolinePlots(ranksPerLength, lengths, Titles.VIOLINE_RANKS_PLOT_2, 
                                   TRUE, TRUE, Extensions.PDF, palette, Strings.LENGTH);
}


#########################################################
#' Reads in the scores i.e. distances and computes the ranks.
#' The scores are then sorted in ascending order 
#' with the help of the sorted ranks.
#' Then for different score-thresholds d in [d_min, d_max]
#' it is looked how many ranks are in the top t ranks
#' or not there
#' under the given threshold d i.e. 
#' every score above the thresold is cut away.
#' The two numbers of how many there (recall) or not there (fallout)
#' are then plotted in percentage into the same plot,
#' where as for the recall the y-axis is used
#' and for the fallout the x-axis.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @export
Experiment1.createROC <- function(passPathGenericName, scoresFolderGenericName, passes, scoreTypes) {
 
  combinedDataPerScoreType <- Encoder.createRocPlotData(passPathGenericName, scoresFolderGenericName, 
                                                        passes, scoreTypes); 
  # plot
  data <- Plotter.plotROC(combinedDataPerScoreType, scoreTypes, Titles.CHARACTERISTIC_5, 
                          TRUE, TRUE, Extensions.PNG);
}


#########################################################
#' Creates a back-to-back histogram for a number of worst and best ranks.
#' 
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param rankIndicesToLookAt {vector} the indices to look at from positive and negative site
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @export
Experiment1.createBackToBackHistogram <- function(passPathGenericName, scoresFolderGenericName, 
                                                  rankIndicesToLookAt, passes) {
  
  samplesData <- Encoder.getRanksAndFilenamesPerPass(passPathGenericName, scoresFolderGenericName, passes);
  
  ranksPerPass <- samplesData$ranksPerPass;
  fileNamesPerPass <- samplesData$fileNamesPerPass;
  
  # plot
  Plotter.plotBackToBackHistogram(ranksPerPass, fileNamesPerPass, 
                                  unlist(rankIndicesToLookAt), TRUE,
                                  TRUE, FALSE);
}