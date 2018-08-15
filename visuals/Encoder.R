#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

ENCODER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Contains methods to encode data for visuals like plotters.

#########################################################
#' Returns the scores for the pattern and its subpatterns at the correct position.
#'
#' @param scoreListsPerPatternLength {list} the scores per pattern length over all possible positions
#' @param position {numerical} the position from which the score should be extracted
#' @param normalize {boolean} tells if scores should be normalized or not (divided by its number of profiles)
#'
#' @return {list} the scores for the given position
#' @export
Encoder.getCorrectPatternPositionScores <- function(scoreListsPerPatternLength, position, normalize) {
  patternLenghts <- length(scoreListsPerPatternLength);
  
  scores <- list();
  
  # iterate over all pattern lengths
  for (i in 1:patternLenghts) {
    score <- scoreListsPerPatternLength[[i]][[position]];
    
    if (normalize) {
      score <- score / i;
    }
    
    scores <- rlist::list.append(scores, score);
  }
  
  return(scores);
}


#########################################################
#' Creates a violine plot data out of the ranks from several computation passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param decreasing {logical} tells if the highest value should have lowest rank
#' @param exception {logical} tells if in the first round the decreasing parameter should be FALSE
#'
#' @return {list} list of ranks per score type from all passes
#' @export
Encoder.createViolinePlotData <- function(passPathGenericName, scoresFolderGenericName, 
                                          passes, scoreTypes, decreasing,
                                          exception) {
  # (C): compute ranks out of scores
  ranksPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    ranksPerScoreType <- list();
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # iterate over each score type
    for (j in 1:length(scoreTypes)) {
      # get score-type
      scoreType <- scoreTypes[[j]];
      
      # read in score-type scores per sample
      scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                   scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
      
      # get scores for correct position of the samples
      scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                               masterChronology$years);
      
      # for each sample: compute the rank of the correct-position score within all scores of the sample
      if (exception && j == 1) {
        ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, FALSE);
      } else {
        ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, decreasing);
      }
      
      # store the ranks per score-type
      ranksPerScoreType <- rlist::list.append(ranksPerScoreType, ranks);
    }
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerScoreType);
  }
  
  combinedRanksPerScoreType <- Encoder.createCombinedRanks(ranksPerPass, passes, scoreTypes);
  
  return(combinedRanksPerScoreType);
}


#########################################################
#' Combines the the ranks from several passes per score type.
#'
#' @param ranksPerPass {list} the list with lists which are showing ranks per score type
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' 
#' @return {list} the ranks of multiple passes added together to a list of ranks per score type
Encoder.createCombinedRanks <- function(ranksPerPass, passes, scoreTypes) {
  ranksPerScoreType <- list();
  
  # iterate over all score types
  for (i in 1:length(scoreTypes)) {
    
    allRanksForScoreType <- list();
    
    # iterate over all passes
    for(j in 1:length(passes)) {
      ranks <- ranksPerPass[[j]][[i]];
      allRanksForScoreType <- rlist::list.append(allRanksForScoreType, unlist(ranks));
    }
    
    ranksPerScoreType <- rlist::list.append(ranksPerScoreType, unlist(allRanksForScoreType, recursive = FALSE));
  }
  
  return(ranksPerScoreType);
}


#########################################################
#' Creates series-length violine plot data out of the ranks from several computation passes.
#'
#' @param seriesLengthPaths {list} the paths to the different series lengths
#' @param lengths {list} the lengths of that paths
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' Hint: It can be be specified more than available.
#' @param scoreType {string} the score type at which to look at
#' @param palette {vector} the colors palette which should be used
#' 
#' @return {list} the ranks of multiple passes added together to a list of ranks per score type
#' @export
Encoder.createSeriesLengthViolinePlotData <- function(seriesLengthPaths, lengths, 
                                                      passPathGenericName, passes, 
                                                      scoreType, palette) {
  # compute ranks out of scores
  ranksPerLength <- list();
  
  # iterate over all series lengths
  for (i in 1:length(seriesLengthPaths)) {
    
    ranksPerPass <- list();
    # iterate over all passes for this length
    for (j in 1:length(passes)) {  # hint: correspondence is not lost -> separate files loaded in correct order (tested)
      passesPath <- paste(seriesLengthPaths[[i]], passes[[j]], sep = Symbols.EMPTY);
      existentPath <- dir.exists(passesPath);
      
      if (existentPath) {
        sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH),
                                          Symbols.REAL_DATA_SEPARATOR);
        masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY),
                                                  Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
        
        scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, Files.SCORES_FOLDER_NAME,
                                                     scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
        
        # get scores for correct position of the samples
        scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                                 masterChronology$years);
        
        # for each sample: compute the rank of the correct-position score within all scores of the sample
        ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, decreasing = FALSE);
        
        # store the ranks per score-type
        ranksPerPass <- rlist::list.append(ranksPerPass, ranks);
      }
    }
    
    ranksPerLength <- rlist::list.append(ranksPerLength, ranksPerPass);
  }
  
  return(ranksPerLength);
}


#########################################################
#' Returns data for different plotTypes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param pass {numerical} the pass which has to be used for computation
#' @param ranksToLookAt {list} the ranks for which the Venn Diagrams should be created 
#' (negative ranks indices gives you worst ranks i.e. -1 = last rank)
#' @param lookOnIndices {logical} tells if to look on indices or ranks i.e. changes ranksToLookAt
#' 
#' @return {list(ranksPerScoreType, filenames)} the ranks per score type and the corresponding sample names
#' @export
Encoder.createVennDiagramData <- function(passPathGenericName, scoresFolderGenericName, 
                                          pass, ranksToLookAt, lookOnIndices) {
  passPath <- paste(passPathGenericName, pass, sep = Symbols.EMPTY);
  
  sampleData <- Loader.readInCurves(paste(passPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH),
                                    Symbols.REAL_DATA_SEPARATOR);
  
  masterChronology <- Loader.readInProfiles(paste(passPath, Paths.PATH ,sep = Symbols.EMPTY),
                                            Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
  
  return(Encoder.__getRanksPerScoreTypeAndSampleNames(scoresFolderGenericName, passPath, 
                                                      sampleData, masterChronology));
}


#########################################################
#' Returns data for different plotTypes.
#'
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passPath {string} the path to the pass
#' @param sampleData {list} the curves for which ranks are computed
#' @param masterChronology {list} the master chronology from which the years are extracted
#' 
#' @return {list(ranksPerScoreType, filenames)} the ranks per score type and the corresponding sample names
Encoder.__getRanksPerScoreTypeAndSampleNames <- function(scoresFolderGenericName, passPath, 
                                                         sampleData, masterChronology) {
  ranksPerScoreType <- list();
  scoresPerSample <- list();
  
  # iterate over all score types
  for (i in 1:length(Parameters.SCORE_TYPES_TO_LOOK_AT)) {
    # get score-type
    scoreType <- Parameters.SCORE_TYPES_TO_LOOK_AT[[i]];
    
    # read in score-type scores per sample
    scoresPerSample <- Loader.readInScores(paste(passPath, Paths.PATH, scoresFolderGenericName,
                                                 scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
    
    # get scores for correct position of the samples
    scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample$scoresData, sampleData,
                                                                             masterChronology$years);
    # for each sample: compute the rank of the correct-position score within all scores of the sample
    ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample$scoresData);
    
    # store the ranks per score-type
    ranksPerScoreType <- rlist::list.append(ranksPerScoreType, ranks);
  }
  
  return(list(ranksPerScoreType = ranksPerScoreType, fileNames = scoresPerSample$fileNames));
}


#########################################################
#' Creates data for an ROC-plot.
#' 
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' 
#' @return {list} the data for an ROC-plot
#' @export
Encoder.createRocPlotData <- function(passPathGenericName, scoresFolderGenericName, passes, scoreTypes) {
  # compute ranks out of scores
  ranksPerPass <- list();
  scoresPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    ranksPerScoreType <- list();
    scoresPerScoreType <- list();
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # iterate over each score type
    for (j in 1:length(scoreTypes)) {
      # get score-type
      scoreType <- scoreTypes[[j]];
      
      # read in score-type scores per sample
      scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                   scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
      
      # get scores for correct position of the samples
      scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                               masterChronology$years);
      
      # for each sample: compute the rank of the correct-position score within all scores of the sample
      ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample);
      
      # store the ranks per score-type
      ranksPerScoreType <- rlist::list.append(ranksPerScoreType, ranks);
      scoresPerScoreType <- rlist::list.append(scoresPerScoreType, scoresForCorrectPositions);
    }
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerScoreType);
    scoresPerPass <- rlist::list.append(scoresPerPass, scoresPerScoreType);
  }
  
  combinedDataPerScoreType <- Encoder.createCombinedData(ranksPerPass, scoresPerPass, passes, scoreTypes);
  
  return(combinedDataPerScoreType);
}


#########################################################
#' Combines the the data from several passes per score type.
#'
#' @param ranksPerPass {list} the list with lists which are showing ranks per score type
#' @param scoresPerPass {list} the list with lists which are showing scores per score type
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' 
#' @return {list(ranksPerScoreType, scoresPerScoreType)} the data of multiple passes added together to a list of ranks per score type
#' @export
Encoder.createCombinedData <- function(ranksPerPass, scoresPerPass, passes, scoreTypes) { 
  ranksPerScoreType <- list();
  scoresPerScoreType <- list();
  
  # iterate over all score types
  for (i in 1:length(scoreTypes)) {
    
    allRanksForScoreType <- list();
    allScoresForScoreType <- list();
    
    # iterate over all passes
    for(j in 1:length(passes)) {
      ranks <- ranksPerPass[[j]][[i]];
      allRanksForScoreType <- rlist::list.append(allRanksForScoreType, unlist(ranks));
      
      scores <- scoresPerPass[[j]][[i]];
      allScoresForScoreType <- rlist::list.append(allScoresForScoreType, unlist(scores));
    }
    
    ranksPerScoreType <- rlist::list.append(ranksPerScoreType, unlist(allRanksForScoreType, recursive = FALSE));
    scoresPerScoreType <- rlist::list.append(scoresPerScoreType, unlist(allScoresForScoreType, recursive = FALSE));
  }
  
  return(list(ranksPerScoreType = ranksPerScoreType, 
              scoresPerScoreType = scoresPerScoreType));
}


#########################################################
#' Returns ranks and filenames of samples per pass.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#'
#' @return {list(ranksPerPass, fileNamesPerPass)} ranks and filenames per pass
#' @export
Encoder.getRanksAndFilenamesPerPass <- function(passPathGenericName, scoresFolderGenericName, passes) {
  ranksPerPass <- list();
  fileNamesPerPass <- list();
  
  # iterate over all passes
  for (i in 1:length(passes)) {
    # retrieve data per pass
    pass <- passes[[i]];
    passPath <- paste(passPathGenericName, pass, sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH),
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passPath, Paths.PATH ,sep = Symbols.EMPTY),
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # get ranks and filenames
    data <- Encoder.__getRanksPerScoreTypeAndSampleNames(scoresFolderGenericName, passPath, 
                                                         sampleData, masterChronology);
    
    ranksPerScoreType <- data$ranksPerScoreType;
    fileNamesOfPass <- data$fileNames;
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerScoreType);
    fileNamesPerPass <- rlist::list.append(fileNamesPerPass, fileNamesOfPass);
  }
  
  return(list(ranksPerPass = ranksPerPass, fileNamesPerPass = fileNamesPerPass));
}


#########################################################
#' Creates ratio bar plot data.
#' 
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param rankIndicesToLookAt {vector} the indices to look at from positive and negative site
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param startYear {number} the year at which the histogram should start
#' @param endYear {number} the year at which the histogram should end
#' 
#' @return {list} the data for a ratio bar-plot
#' @export
Encoder.createRatioBarPlotData <- function(passPathGenericName, scoresFolderGenericName, 
                                           rankIndicesToLookAt, passes, startYear, endYear) {
  # retrieve data
  samplesData <- Encoder.getRanksAndFilenamesPerPass(passPathGenericName, scoresFolderGenericName, passes);
  
  ranksPerPass <- samplesData$ranksPerPass;
  fileNamesPerPass <- samplesData$fileNamesPerPass;
  
  # create data
  positiveYears <- Encoder.getYearHistogramEncoding(ranksPerPass, fileNamesPerPass,
                                                    rankIndicesToLookAt, TRUE);
  
  negativeYears <- Encoder.getYearHistogramEncoding(ranksPerPass, fileNamesPerPass, 
                                                    -rankIndicesToLookAt, TRUE);
  
  # encode
  ratios <- list();  # ratios 
  
  for (i in startYear:endYear) {
    numberOfPositiveYears <- length(positiveYears[positiveYears == i]) + 1;  
    numberOfNegativeYears <- length(negativeYears[negativeYears == i]) + 1;  # +1 to avoid division by zero
    
    ratio <- numberOfPositiveYears / numberOfNegativeYears;
    if (ratio <= 1) ratio = 0;
    
    ratios <- rlist::list.append(ratios, ratio);
  }
  
  return(ratios);
}


#########################################################
#' Creates ratio bar plot data.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param palette {vector} the colors palette which should be used
#' @param shapes {vector} the shapes palette which should be used
#' @param weighted {logical} tells if the lowest score should be weighted with its delta score or not  
#'
#' @return {list} creates the data for a scatter difference plot
#' @export
Encoder.createDifferenceScatterPlotData <- function(passPathGenericName, scoresFolderGenericName, 
                                                    passes, scoreType, palette, shapes, weighted) {
  ranksPerPass <- list();
  # correctScoresPerPass <- list();
  rank1ScorePerPass <- list();
  rank2ScorePerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # read in score-type scores per sample
    scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                 scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
    
    # get prediction scores (x-axis)
    rank1PredictionScoresPerSample <- Analyzer.getRankScore(scoresPerSample, 1);
    rank2PredictionScoresPerSample <- Analyzer.getRankScore(scoresPerSample, 2);
    
    if (weighted) {  # change scoresPerSample rank 1 score depending on distance between rank 1 and rank 2 score
      # compute weights
      weights <- unlist(rank2PredictionScoresPerSample) - unlist(rank1PredictionScoresPerSample);
      
      scoresPerSample <- Analyzer.changeSampleScore(scoresPerSample, rank1PredictionScoresPerSample, weights);
      
      # rank 1/rank 2 prediction scores has to be recomputed
      rank1PredictionScoresPerSample <- Analyzer.getRankScore(scoresPerSample, 1);
      rank2PredictionScoresPerSample <- Analyzer.getRankScore(scoresPerSample, 2);
    }
    
    # get scores for correct position of the samples
    scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                             masterChronology$years);
    
    # for each sample: compute the rank of the correct-position score within all scores of the sample
    ranksPerSample <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample);
    
    # store the ranks per score-type
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerSample);
    
    # correctScoresPerPass <- rlist::list.append(correctScoresPerPass, scoresForCorrectPositions);
    rank1ScorePerPass <- rlist::list.append(rank1ScorePerPass, rank1PredictionScoresPerSample);
    rank2ScorePerPass <- rlist::list.append(rank2ScorePerPass, rank2PredictionScoresPerSample);
  }
  
  # encoding
  xData <- list();
  yData <- list();
  group <- list();
  
  # iterate over each pass
  for (i in 1:length(rank1ScorePerPass)) {  
    rank1PredictionScoresPerSample <- rank1ScorePerPass[[i]];
    rank2PredictionScoresPerSample <- rank2ScorePerPass[[i]];
    ranksPerSample <- ranksPerPass[[i]];
    
    # iterate over each prediction
    for (j in 1:length(rank1PredictionScoresPerSample)) {
      predRank1Score <- rank1PredictionScoresPerSample[[j]];
      predRank2Score <- rank2PredictionScoresPerSample[[j]];
      rank <- ranksPerSample[[j]];
      
      xScore <- predRank2Score - predRank1Score;
      
      xData <- rlist::list.append(xData, xScore);
      yData <- rlist::list.append(yData, rank);
      group <- rlist::list.append(group, Strings.POINT);
    }
  }
  
  return(list(xData = xData, yData = yData, group = group));
}


#########################################################
#' Creates the box plots data for the histogram based voting approach
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
#'
#' @return {list} the data for the box plot
#' @export
Encoder.computeBoxDataPlot <- function(passPathGenericName, scoresFolderGenericName,
                                       passes, scoreTypes, palette, topYearsCount, 
                                       powersetApproach, doubleWeighting,
                                       logarithmicWeighting, minimumLength) {
  # compute ranks out of scores
  ranksPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    ranksPerScoreType <- list();
    
    # iterate over each score type
    for (j in 1:length(scoreTypes)) {
      # get score-type
      scoreType <- scoreTypes[[j]];
      
      # read in score-type scores per sample
      scoresPerSample <- Loader.readInIndividualScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                             scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
      
      topYearsPerSample <- Analyzer.getTopYearsPerSample(scoresPerSample$scoresTables, topYearsCount, 
                                                         powersetApproach, doubleWeighting, logarithmicWeighting, 
                                                         FALSE, minimumLength);
      
      ranks <- Analyzer.getPeakRanks(sampleData, topYearsPerSample);
      
      ranksPerScoreType <- rlist::list.append(ranksPerScoreType, ranks);
    }
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerScoreType);
  }
  
  combinedRanksPerScoreType <- Encoder.createCombinedRanks(ranksPerPass, passes, scoreTypes);
  
  return(combinedRanksPerScoreType);
}


#########################################################
#' Creates the box plots data for the histogram based voting approach
#' out of the ranks from several computation passes.
#'
#' @param number {numerical} the number of samples on which you want to look on
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the score type name of the ring-width based method
#' 
#' @return {list} the data for the rank enhancement plot
#' @export
Encoder.createRankEnhancementPlotData <- function(number, passPathGenericName, 
                                                  scoresFolderGenericName, 
                                                  passes, scoreType) {
  
  # get top rated samples from ring-width based approach
  topSamplesAndRanks <- Analyzer.getTopSamples(number, passPathGenericName, 
                                               scoresFolderGenericName, 
                                               passes, scoreType);
  
  widthBasedRanks <- topSamplesAndRanks$ranks;
  
  # look up scores of the same samples under bucket based best method d)
  fileNames <- topSamplesAndRanks$fileNames;
  
  # all scoring-files are encoded in the same way, replace in fullnames the strings
  # such that for the best method d) from bucket-based approach the ranks can be looked up
  fileNames <- gsub(Strings.TO_REPLACE_3, Strings.REPLACEMENT_3, fileNames);
  fileNames <- gsub(Strings.TO_REPLACE_4, Strings.REPLACEMENT_4, fileNames);
  fileNames <- gsub(Strings.TO_REPLACE_5, Strings.REPLACEMENT_5, fileNames);
  fileNames <- gsub(Strings.TO_REPLACE_6, Strings.REPLACEMENT_3, fileNames);
  
  # look up the ranks of the files above under the bucket-based approach best method
  bucketBasedRanks <- Analyzer.getFilesRanks(fileNames);
  
  # encode plot data
  x <- 1:length(bucketBasedRanks);
  y <- list();
  g <- list();
  
  y <- rlist::list.append(y, widthBasedRanks);
  g <- rlist::list.append(g, rep(Strings.RING_WIDTH_BASED, length(widthBasedRanks)));
  y <- rlist::list.append(y, bucketBasedRanks);
  g <- rlist::list.append(g, rep(Strings.BUCKET_BASED, length(bucketBasedRanks)));
  
  return(list(x = x, y = y, g = g));
}


#########################################################
#' Creates the rank histogram data after applying the buckets-based
#' approach on the top years from the ring-width based approach.
#'
#' @param topYears {numerical} the number of top scores that have to be retrieved 
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType of the ring-width approach
#' @param histogramScoreType {list} the scoreType of the buckets-based method
#' 
#' @return {vector} the data for a rank histogram plot
#' @export
Encoder.createsRankHistogramForTopYears <- function(topYears, passPathGenericName, 
                                                    scoresFolderGenericName, passes, 
                                                    scoreType, histogramScoreType) {
  
  bucketBasedRanks <- Analyzer.getCorrespondingRanksAndOutcasts(topYears, passPathGenericName, 
                                                                scoresFolderGenericName, passes, 
                                                                scoreType, histogramScoreType);
  yData <- unlist(bucketBasedRanks$ranks);
  
  print(Strings.OUTLIERS);
  print(bucketBasedRanks$outcasts);
  
  return(yData);
}


#########################################################
#' Creates data for profiles per year histogram.
#'
#' @param startYear {number} the year at which the histogram should start
#' @param endYear {number} the year at which the histogram should end
#' 
#' @return {vector} the years
#' @export
Encoder.createProfilesPerYearData <- function(startYear, endYear) {
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  perYearProfiles <- CurvesMiner.collectData(curvesData)$profilesPerYear;
  
  years <- list();  # stores all years
  
  # encode for histogram
  for (i in startYear:endYear) {
    numberOfProfiles <- length(perYearProfiles[[toString(i)]]);  # how often to insert this year
    years <- rlist::list.append(years, rep(i, numberOfProfiles));
  }
  
  return(years);
}


#########################################################
#' Creates an all-vs-all plot data for the distances between profiles.
#'
#' @param startYear {number} the year at which to start
#' @param endYear {number} the year at which to end
#' @param meanDistances {logical} tells if the mean-distances should be plotted
#' @param title {string} the title for the plot
#' 
#' @return {list} data for a all-vs-all plot
#' @export
Encoder.createAllVsAllDistancesPlotData <- function(startYear, endYear, meanDistances = FALSE, 
                                                    title = Titles.PROFILES_DISTANCES_PER_YEAR) {
  distancesList <- list();
  
  # iterate over all years
  for (i in startYear:endYear) {
    path <- paste(Paths.ALL_VS_ALL_DISTANCES, i, Extensions.CSV, sep = Symbols.EMPTY);
    mat <- Loader.readInMatrix(path, Symbols.REAL_DATA_SEPARATOR);
    
    if (meanDistances) {
      numMeanDistances <- dim(mat)[1] - 1;  # n profiles -> n-1 distances
      distances <- colSums(mat) / numMeanDistances;
    } else {
      distances <- mat[lower.tri(mat)];
    }
    
    distancesList <- rlist::list.append(distancesList, distances);
  }
  
  return(distancesList);
}


#########################################################
#' Creates a Venn diagram for the ranks of 4 measurement types.
#' So it is looked if the sets of samples for some rank have an overlap.
#'
#' @param ranksPerScoreType {list} list of sublists, 
#' where every sublist contains the ranks of the different samples
#' @param fileNames {list} the names of the files from the ranks
#' (maximum number of different score types = length(Palettes.NATURE))
#' @param ranksToLookAt {vector} the ranks for which the plot should be created
#' @param pass {numerical} the pass to encode into file names
#' @param save {logical} save the with a list  plot or not
#' 
#' @return {list} the data for a Venn diagram
#' @export
Encoder.creatVennDiagram1Data <- function(ranksPerScoreType, fileNames, ranksToLookAt, pass, save) {
  fileNamesPerSet <- list();
  
  # iterate over each score type
  for (i in 1:length(ranksPerScoreType)) {
    ranksOfScoreType <- ranksPerScoreType[[i]];
    # get all indices from ranks at which you want look at
    indices <- which(unlist(ranksOfScoreType) %in% ranksToLookAt);
    names <- fileNames[indices];
    fileNamesPerSet <- rlist::list.append(fileNamesPerSet, unlist(names));
  }
  
  return(fileNamesPerSet);
}


#########################################################
#' Creates Venn diagram data for the ranks of 4 measurement types.
#' So it is looked if the sets of samples for some rank have an overlap.
#'
#' @param ranksPerScoreType {list} list of sublists, 
#' where every sublist contains the ranks of the different samples
#' @param fileNames {list} the names of the files from the ranks
#' (maximum number of different score types = length(Palettes.NATURE))
#' @param rankIndicesToLookAt {vector} the rank indices in ascending ordered ranks for which the Venn Diagrams should be created 
#' (negative ranks indices gives you worst ranks i.e. -1 = last rank)
#' @param print {logical} show the plot or not
#' @param pass {numerical} the pass to encode into file names
#' @param save {logical} save plot and a list with maximum intersection elements or not
#' 
#' @return {list} the data for a Venn diagram
#' @export
Encoder.creatVennDiagram2Data <- function(ranksPerScoreType, fileNames, rankIndicesToLookAt, pass, save) {
  fileNamesPerSet <- list();
  
  # iterate over each score type
  for (i in 1:length(ranksPerScoreType)) {
    ranksOfScoreType <- ranksPerScoreType[[i]];
    
    # get sorted filenames
    sortedRanksIndices <- order(unlist(ranksOfScoreType));
    sortedNames <- fileNames[sortedRanksIndices];
    
    # convert negative rank indices into positive
    rankIndicesToLookAt[rankIndicesToLookAt < 0] <- rankIndicesToLookAt[rankIndicesToLookAt < 0] + length(sortedNames) + 1;  
    
    # get names of corresponding rank indices in sorted names (+1 because else last element can't be selected)
    names <- sortedNames[rankIndicesToLookAt];
    
    fileNamesPerSet <- rlist::list.append(fileNamesPerSet, unlist(names));
  }
  
  return(fileNamesPerSet);
}


#########################################################
#' Returns the histogram encoding for a year histogram for given ranks.
#'
#' @param ranksPerScorePass {list} list of ranks per score type
#' @param fileNamesPerPass {list} the names fo the files per pass
#' @param rankIndicesToLookAt {vector} the rank indices in ascending ordered ranks for which the Venn Diagrams should be created 
#' (negative ranks indices gives you worst ranks i.e. -1 = last rank)
#' @param maxIntersection {logical} if you want create a plot for the intersection of all four methods
#'
#' @return {yearsList} the histogram years
#' @export
Encoder.getYearHistogramEncoding <- function(ranksPerScorePass, fileNamesPerPass, 
                                             rankIndicesToLookAt, maxIntersection) {
  # what to plot
  yearsList <- list();
  
  # encoding
  for (j in 1:length(ranksPerScorePass)) {  # iterate over all passes
    ranksPerScoreType <- ranksPerScorePass[[j]];
    files <- fileNamesPerPass[[j]];
    
    rankRoundIndicesToLookAt <- rankIndicesToLookAt;  # create copy
    
    # iterate over each score type    
    fileNamesPerSet <- list(); 
    
    # print("names:");
    for (i in 1:length(ranksPerScoreType)) {
      ranksOfScoreType <- ranksPerScoreType[[i]];
      
      # get sorted filenames with the help of the ranks
      sortedRanksIndices <- order(unlist(ranksOfScoreType));
      sortedNames <- files[sortedRanksIndices];
      
      # convert negative rank indices into positive
      rankRoundIndicesToLookAt[rankRoundIndicesToLookAt < 0] <- rankRoundIndicesToLookAt[rankRoundIndicesToLookAt < 0] + length(sortedNames) + 1;  
      
      # get filenames or corresponding rank indices in sorted names (+1 because else last element can't be selected)
      names <- sortedNames[rankRoundIndicesToLookAt];
      # print(length(names));
      fileNamesPerSet <- rlist::list.append(fileNamesPerSet, unlist(names));
    }
    
    # get the files
    fileNames <- Symbols.EMPTY;
    
    # print("individual:");
    # print(length(fileNamesPerSet[[1]]));
    # print(length(fileNamesPerSet[[2]]));
    # print(length(fileNamesPerSet[[3]]));
    # print(length(fileNamesPerSet[[4]]));
    
    if (maxIntersection) {  # data from intersection of all filenames
      fileNames <- Reduce(intersect, fileNamesPerSet[1:4]);
    } else {
      fileNames <- unlist(fileNamesPerSet);
    }
    
    # get the years
    fileNames <- gsub(Strings.TO_REPLACE_1, Symbols.EMPTY, as.list(fileNames));
    yearsStrings <- gsub(Strings.TO_REPLACE_2, Strings.REPLACEMENT, as.list(fileNames));
    
    print(length(yearsStrings));
    
    if (length(yearsStrings) > 0) {
      # iterate over each string "1935-1944" -> and create sequences 1935:1944
      for (k in 1:length(yearsStrings)) {
        yearString <- stringr::str_extract(yearsStrings[[k]], Expression.YEARS);
        startYear <- stringr::str_sub(yearString, 1, nchar(yearString) - 5);
        endYear <- stringr::str_sub(yearString, 6, nchar(yearString));
        
        years <- startYear:endYear;
        yearsList <- rlist::list.append(yearsList, unlist(years));
      }
    }
  }
  
  return(unlist(yearsList));
} 


#########################################################
#' Creates data for a ROC-plot.
#'
#' @param dataPerScoreType {list(ranksPerScoreType, scoresPerScoreType)} two lists of sublists, 
#' where every sublist contains the ranks (respectively the score) of the different samples
#' @param types {list} list of possible score types 
#' @param maxTopRank {numerical} the maximum number of top interval [1, 2, .., maxTopRank]
#' 
#' @return {list(truePositives, falsePositives, tresholdValues, type)} the data for the ROC-plot
Encoder.createRocPlotEncoding <- function(dataPerScoreType, types, maxTopRank = Parameters.TOP_RANK) {
  # encode data
  ranksPerScoreType <- dataPerScoreType$ranksPerScoreType;
  scoresPerScoreType <- dataPerScoreType$scoresPerScoreType;
  
  falsePositives <- list();  # x-axis
  truePositives <- list();  # y-axis
  tresholdValues <- list();  # marker on curve
  type <- list();
  
  # iterate over each score type
  for (i in 1:length(types)) {  
    ranks <- unlist(ranksPerScoreType[[i]]);
    scores <- unlist(scoresPerScoreType[[i]]);
    
    sortedScores <- sort(scores);
    
    # retrieve (d_S | r_S \in {1,..,t}) and (d_S | r_S \notin {1,..,t})
    relevant <- length(ranks[ranks <= maxTopRank]);
    irrelevant <- length(ranks[ranks > maxTopRank]);
    
    # use every score in sortedScores as a threshold d
    for (j in 1:length(sortedScores)) {
      threshold <- sortedScores[[j]];  # d
      
      # retrieve (d_S | d_S < d) & (d_S | r_S \in {1,..,t}) and ..
      indicesWithScoresLowerThreshold <- which(scores < threshold);
      ranksWithScoresLowerThreshold <- ranks[indicesWithScoresLowerThreshold];
      
      retrievedAndRelevant <- length(ranksWithScoresLowerThreshold[ranksWithScoresLowerThreshold <= maxTopRank]);
      retrievedAndIrrelevant <- length(ranksWithScoresLowerThreshold[ranksWithScoresLowerThreshold > maxTopRank]);
      
      recall <- retrievedAndRelevant/relevant;  # (d_S | d_S < d) & (d_S | r_S \in {1,..,t})
      fallout <- retrievedAndIrrelevant/irrelevant;  # (d_S | d_S < d) & (d_S | r_S \notin {1,..,t})
      
      falsePositives <- rlist::list.append(falsePositives, fallout);
      truePositives <- rlist::list.append(truePositives, recall);
      tresholdValues <- rlist::list.append(tresholdValues, round(threshold, 2));
      type <- rlist::list.append(type, i);
    }
  }
  
  return(list(truePositives = truePositives, falsePositives = falsePositives, 
              tresholdValues = tresholdValues, type = type));
}


#########################################################
#' Creates the encoding for the correlation coefficient histogram.
#'
#' @param consensiPath {string} the path to the consensi
#' @param consensusName {string} the name of the consensus
#'
#' @return {list(xData, yData)} the x- and y-data vectors for the histogram
#' @export
Encoder.createCorrelationCoefficientHistogramEncoding <- function(consensiPath, consensusName) {
  masterChronologyData <- Loader.readInProfiles(consensiPath, consensusName, Symbols.REAL_DATA_SEPARATOR);
  
  # encode x data
  xData <- list();
  
  for (i in 1:(length(masterChronologyData$years)-1)) {  # -1, since between n-values are n-1 differences
    value <- as.numeric(masterChronologyData$years[[i]]) + 0.5;  # +0.5 since between two years
    xData <- rlist::list.append(xData, value);
  }
  
  xData <- unlist(xData);  
  neighbourScores <- Analyzer.compareNeighbouredProfiles(masterChronologyData$profiles);
  
  return(list(xData = xData, yData = unlist(neighbourScores)));
}


#########################################################
#' Creates the encoding for the correlation coefficient histogram.
#'
#' @param consensiPath {string} the path to the consensi
#' @param consensusName {string} the name of the consensus
#' @param correlationCoefficientLimit {numerical} the limit below which all years defined as pointer-years
#' for each sample with the ring-width based approach
#' @param correlationCoefficientLimitForEvents {numerical} the limit below which years defined as event-years
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType of the buckets-based method
#'
#' @return {list(yData, outRanks, nonCharacteristicRanks, characteristicRanks)} the encoding data
#' @export
Encoder.createRanksPointerHeuristicHistogramEncoding <- function(consensiPath, consensusName,
                                                                 correlationCoefficientLimit,
                                                                 correlationCoefficientLimitForEvents,
                                                                 passPathGenericName, 
                                                                 scoresFolderGenericName, passes, 
                                                                 scoreType) {
  # load
  masterChronologyData <- Loader.readInProfiles(consensiPath, consensusName, Symbols.REAL_DATA_SEPARATOR);
  
  # determine pointer-years
  allYears <- unlist(masterChronologyData$years);
  neighbourScores <- Analyzer.compareNeighbouredProfiles(masterChronologyData$profiles);
  
  # starting with value 2, since in the first year pointer-years are non-detectable
  pointerYears <- Analyzer.getCharacteristicYears(neighbourScores,
                                                  unlist(allYears[2:length(allYears)]),
                                                  correlationCoefficientLimit);
  
  allRanks <- list();
  numOutCasts <- 0;
  nonCharacteristicRanks <- list();
  characteristicRanks <- list();

  # determine samples with event-years
  # iterate over all passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);

    # look which samples contain pointer-years and which not
    detection <- Analyzer.splitIntoCharacteristicYears(sampleData, correlationCoefficientLimitForEvents);
    
    samplesWithoutEventYears <- detection$samplesWithoutCharacteristicYears;
    samplesWithEventYears <- detection$samplesWithCharacteristicYears;
    eventYearIndicesPerSample <- detection$characteristicYearsIndicesPerSample;
    
    scoresPath <- paste(passesPath, Symbols.SLASH, scoresFolderGenericName, 
                        scoreType, Symbols.SLASH, sep = Symbols.EMPTY);
    
    eventFreeRanks <- Analyzer.getRanksForUsualSamples(allYears, scoresPath, samplesWithoutEventYears, scoreType);  # as before
    eventRanksData <- Analyzer.getRanksForCharacteristicSamples(allYears, scoresPath, samplesWithEventYears, eventYearIndicesPerSample,
                                                                scoreType, pointerYears);
    
    eventRanks <- eventRanksData$predictedRanks;
    outcasts <- eventRanksData$outcasts;
    
    # append
    nonCharacteristicRanks <- rlist::list.append(nonCharacteristicRanks, eventFreeRanks);
    characteristicRanks <- rlist::list.append(characteristicRanks, eventRanks);

    combinedRanks <- c(eventFreeRanks, eventRanks);
    allRanks <- rlist::list.append(allRanks, combinedRanks);
    numOutCasts <- numOutCasts + outcasts;
  }
  
  return(list(yData = unlist(allRanks), outcasts = numOutCasts, 
              nonCharacteristicRanks = unlist(nonCharacteristicRanks), 
              characteristicRanks = unlist(characteristicRanks)));
}


#########################################################
#' Creates the encoding for the goodness histogram for delta-peaks.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} contains the scoreType
#' @param topYearsCount {numerical} how many years considered for each profile
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
#' 
#' @return {list(xData, yData, yOutcast)} the encoded data for a plot
#' @export
Encoder.createGoodnessHistogramEncoding <- function(passPathGenericName, scoresFolderGenericName,
                                                    passes, scoreType, topYearsCount, 
                                                    minimumLength) {
  # compute ranks out of scores
  ranksPerPass <- list();
  deltaPeaksPerPass  <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    ranksPerScoreType <- list();
    
    # read in score-type scores per sample
    scoresPerSample <- Loader.readInIndividualScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                           scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
    
    topYearsPerSample <- Analyzer.getTopYearsPerSample(scoresPerSample$scoresTables, topYearsCount, 
                                                       TRUE, FALSE, FALSE, 
                                                       FALSE, minimumLength);
    
    ranksPerSample <- Analyzer.getPeakRanks(sampleData, topYearsPerSample);
    deltaPeaksPerSample <- Analyzer.getWeightedDeltaPeaks(topYearsPerSample);  # quality per sample
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerSample);
    deltaPeaksPerPass <- rlist::list.append(deltaPeaksPerPass, deltaPeaksPerSample);
  }
  
  ranksPerSample <- unlist(ranksPerPass);
  deltaPeaksPerSample <- unlist(deltaPeaksPerPass);
  
  # encoding
  xData <- list();  # quality
  yData <- list();  # correctness
  
  yOutcast <- list();
  
  # iterate over all samples
  for (i in 1:length(ranksPerSample)) {
    rank <- ranksPerSample[[i]];
    delta <- deltaPeaksPerSample[[i]];
    
    if (delta == Defaults.ALL_VOTED_FOR_THE_SAME) {  # if outcast i.e. all voted for the same
      if (rank == 1) {
        yOutcast <- rlist::list.append(yOutcast, 1);
      } else {
        yOutcast <- rlist::list.append(yOutcast, -1);
      }
    } else {  # more than one bar in the histogram
      xData <- rlist::list.append(xData, delta);
      
      if (rank == 1) {
        yData <- rlist::list.append(yData, 1);
      } else {
        yData <- rlist::list.append(yData, -1);
      }
    }
  }
  
  return(list(xData = unlist(xData), yData = unlist(yData), yOutcast = unlist(yOutcast)));
}


#########################################################
#' Creates the encoding for the goodness histogram for delta-scores.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param years {list} the years for which the diagram has to be created
#'
#' @return {list(xData, yData)} the encoded data for a plot
#' @export
Encoder.createGoodnessHistogramEncoding2 <- function(passPathGenericName, 
                                                     scoresFolderGenericName,
                                                     passes, scoreType, years) {
  # compute ranks out of scores
  ranksPerPass <- list();
  deltaScoresPerPass  <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);

    # read in score-type scores per sample
    scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                 scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
    
    # for each sample: compute the rank of the correct-position score within all scores of the sample
    scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                             years);
    # ranks
    ranksPerSample <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample);
    
    rank1PredictionScoresPerSample <- unlist(Analyzer.getRankScore(scoresPerSample, 1));
    rank2PredictionScoresPerSample <- unlist(Analyzer.getRankScore(scoresPerSample, 2));
    
    # deltas
    deltaScoresPerSample <- rank2PredictionScoresPerSample - rank1PredictionScoresPerSample;
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerSample);
    deltaScoresPerPass <- rlist::list.append(deltaScoresPerPass, deltaScoresPerSample);
  }
  
  ranksPerSample <- unlist(ranksPerPass);
  deltaScoresPerPass <- unlist(deltaScoresPerPass);
  
  return(Encoder.__getEncodingData(ranksPerSample, deltaScoresPerPass));
}


#########################################################
#' Returns scores and ranks encoded in x- and y-data.
#'
#' @param scores {vector} the scores which should be encoded
#' @param ranks {vector} the ranks which should be encoded
#'
#' @return {list(xData, yData)} the encoded data for a plot
Encoder.__getEncodingData <- function(ranks, scores) {
  # encoding
  xData <- scores;  # quality
  yData <- list();  # correctness
  
  # iterate over all samples
  for (i in 1:length(ranks)) {
    rank <- ranks[[i]];
    
    if (rank == 1) {
      yData <- rlist::list.append(yData, 1);
    } else {
      yData <- rlist::list.append(yData, -1);
    }
  }
  
  return(list(xData = unlist(xData), yData = unlist(yData)));
}


#########################################################
#' Creates the encoding for a rank/p-value scatter plot.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param years {list} the years for which the diagram has to be created
#' @param threshold {numeric} the threshold at which the p-value has to be computed
#' @param computeDeltaScores {logical} tells if the p-value for delta-scores should be computed (threshold won't work)
#'
#' @return {list(ranksPerSample, pValuesPerSample)} the encoded data for a plot
#' @export
Encoder.createRankPValueEncoding <- function(passPathGenericName, scoresFolderGenericName, 
                                                     passes, scoreTypes, threshold, computeDeltaScores = FALSE) {
  
  data <- Analyzer.getPValuesAndCorrespondingRanks(passPathGenericName, scoresFolderGenericName, 
                                                   passes, scoreTypes, threshold, computeDeltaScores);
  
  # encode
  ranksPerSample <- unlist(data$ranksPerPass);
  pValuesPerSample <- unlist(data$pValuesPerPass);
  
  return(list(ranksPerSample = ranksPerSample, pValuesPerSample = pValuesPerSample));
}



#########################################################
#' Creates the encoding for the goodness histogram for p-values.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param years {list} the years for which the diagram has to be created
#' @param threshold {numeric} the threshold at which the p-value has to be computed
#' @param computeDeltaScores {logical} tells if the p-value for delta-scores should be computed (threshold won't work)
#'
#' @return {list(xData, yData)} the encoded data for a plot
#' @export
Encoder.createGoodnessHistogramEncoding3 <- function(passPathGenericName, scoresFolderGenericName, 
                                                     passes, scoreTypes, threshold, computeDeltaScores = FALSE) {
  
  data <- Analyzer.getPValuesAndCorrespondingRanks(passPathGenericName, scoresFolderGenericName, 
                                                   passes, scoreTypes, threshold, computeDeltaScores);
  
  # encode
  ranksPerSample <- unlist(data$ranksPerPass);
  pValuesPerSample <- unlist(data$pValuesPerPass);
  
  return(Encoder.__getEncodingData(ranksPerSample, pValuesPerSample));
}


#########################################################
#' Creates the encoding for the goodness histogram for p-values.
#'
#' @param startYear {numerical} the year at which the histogram should start
#' @param endYear {numerical} the year at which the histogram should end
#' @param perYearProfiles1 {list} the profiles per Year before some action on them
#' @param perYearProfiles2 {list} the profiles per Year after some action on them
#'
#' @return {list(yearsBefore, yearsAfter)} the years before and after
#' @export
Encoder.createHistogramsEncoding <- function(startYear, endYear, perYearProfiles1, perYearProfiles2) {
  yearsBefore <- list();  # stores all years
  yearsAfter <- list();
  
  # encode for histogram
  for (i in startYear:endYear) {
    numberOfProfilesBefore <- length(perYearProfiles1[[toString(i)]]);  # how often to insert this year
    numberOfProfilesAfter <- length(perYearProfiles2[[toString(i)]]);  # how often to insert this year
    yearsBefore <- rlist::list.append(yearsBefore, rep(i, numberOfProfilesBefore));
    yearsAfter <- rlist::list.append(yearsAfter, rep(i, numberOfProfilesAfter));
  }
  
  return(list(yearsBefore = yearsBefore, yearsAfter = yearsAfter));
}


#########################################################
#' Return the mean encoding for goodness plots like p-value or delta-scores plots
#' against the rank truth mean.
#' 
#' @param xData {vector} the xData from first goodness plot (e.g. delta-score against rank)
#' @param yData {vector} the yData from first goodness plot
#' @param range {numeric} the range in which you want measure the mean
#' @param incrementFrequency {numeric} how often the mean has to be measured
#'
#' @return {list(xData, yData)} the encoded data
#' @export
Encoder.getMeanRangeEncoding <- function(xData, yData, range, incrementFrequency) {
  startValue <- 0;
  nextValue <- range;
  
  xNewData <- list();
  yNewData <- list();
  
  for (i in 1:incrementFrequency) {
    value <- which(xData > startValue & xData <= nextValue);
    average <- mean(yData[value]);
    center <- range/2 + range * (i-1);
    
    # print
    print(length(yData[value]));
    
    # new value
    xNewData <- rlist::list.append(xNewData, center);
    yNewData <- rlist::list.append(yNewData, average);
    
    startValue <- startValue + range;
    nextValue <- nextValue + range;
  }
  
  # encode
  xData <- unlist(xNewData);
  yData <- unlist(yNewData);
  
  return(list(xData = xData, yData = yData));
}


#########################################################
#' Creates a histogram storing the number of trees per given rank.
#' If not specified then a distribution for all trees is computed.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param passPathGenericNameCounts {string} the generic name of the pass path for the counts (without the pass number) 
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} the scoreType for which the histogram should be computed
#' @param rankLookAt {numeric} the rank at which it should be looked at
#'
#' @return {vector} counts of the rank i rated samples
#' @export
Encoder.createHistogramPerSampleTreesEncoding <- function(passPathGenericName, passPathGenericNameCounts,
                                                          scoresFolderGenericName, passes, scoreType, rankLookAt = 0) {
  # counts per pass
  countsPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    passesPathCounts <- paste(passPathGenericNameCounts, passes[[i]], sep = Symbols.EMPTY);
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # read in score-type scores per sample
    scoresPerSample <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                                 scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR)$scoresData;
    
    countsPerSample <- Loader.readInCounts(paste(passesPathCounts, Paths.PATH, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR);
    
    # get scores for correct start years of the samples
    scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                             masterChronology$years);
    
    countsForCorrectPositions <- Analyzer.getCountsForCorrectPositions(countsPerSample, sampleData,  masterChronology$years);
    allCounts <- unlist(countsForCorrectPositions);
    
    # for each sample: compute the rank of the correct-position score within all scores of the sample
    ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, FALSE);
    
    if (rankLookAt != 0) {
      specialRankPositions <- which(ranks == rankLookAt);
      counts <- allCounts[specialRankPositions];
    } else {
      counts <- allCounts;
    }
    
    # store the ranks per score-type
    countsPerPass <- rlist::list.append(countsPerPass, counts);
  }
  
  return(unlist(countsPerPass));
}


#########################################################
#' Extracts the counts of each sample at its correct position in the consensus.
#'
#' @param countsPerSample {list} the scores per sample
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param years {list} the years of the consensus
#'
#' @return {list} score per sample
Analyzer.getCountsForCorrectSamplePositions <- function(countsPerSample, samples, years) {
  counts <- list();
  
  # get for each sample the correct score
  for (i in 1:length(samples)) {
    countsPerSample <- countsPerSample[[i]];
    sample <- samples[[i]];
    correctScore <- Analyzer.__getCorrectPatternPositionScore(scoresOfSample, 
                                                              sample$years[[1]], years);
    counts <- rlist::list.append(counts, correctScore);
  }

  return(counts);
}
