#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

EXPERIMENT_0_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes the test experiments before having any real dataset.

#########################################################
#' Starts the experiment.
#'
#' @param datasetPath {string} the path to the dataset
#' @param qualityName {string} the generic name of the files (without its number)
Experiment0.start <- function(datasetPath, qualityName) {
  # iterate over all possible distance functions
  for (i in 1:Defaults.NUMBER_OF_DISTANCE_FUNCTIONS) {
    distanceFunction <- Experiment0.DISTANCE_FUNCTIONS[[i]];
    
    # iterate over test set
    for (j in 1:Parameters.TEST_SET_SIZE) {
      # get parameters
      treeForConsensus <- Experiment0.CONSENSUS_OF_TREE[[j]];
      subpatternToCutOut <- Experiment0.PATTERNS[[j]];
      patternToLoad <- paste(qualityName, Experiment0.PATTERN_OF_TREE[[j]], sep = Symbols.EMPTY);
      
      #########################################################
      
      # preprocessing: computation of consensus etc. #
      curvesData <- Loader.readInCurves(datasetPath, Symbols.ARTIFICIAL_DATA_SEPARATOR);
      treeData <- CurvesMiner.filterTree(curvesData, treeForConsensus);
      treeConsensus <- CurvesMiner.computeTreeConsensus(treeData, distanceFunction);
      
      #########################################################
      
      # check correctness #
      Experiment0.__correctnessCheck(curvesData, treeData, treeConsensus);
      
      #########################################################
      
      # exercise #
      patternProfilesData <- Loader.readInProfiles(datasetPath, patternToLoad, Symbols.ARTIFICIAL_DATA_SEPARATOR);
      patternProfiles <- CurvesMiner.getSubset(patternProfilesData$profiles, subpatternToCutOut[[1]], subpatternToCutOut[[2]]);
      
      # a)
      scoreListsPerPatternLength <- Analyzer.computeMicaScores(treeConsensus, patternProfiles, distanceFunction);
      
      # b)
      parametersToStore <- list(distanceFunction = distanceFunction,
                                loadedPattern = patternToLoad,
                                subpatternStart = subpatternToCutOut[[1]], 
                                subpatternEnd = (subpatternToCutOut[[1]] + length(scoreListsPerPatternLength)),
                                treeForConsensus = treeForConsensus);
      
      Plotter.plotHistogramsForScores(scoreListsPerPatternLength, subpatternToCutOut[[1]], TRUE, FALSE, parametersToStore);
    }
  }
}


#########################################################
#' Shows the user some information which can be used to see if
#' the data was loaded correctly.
#'
#' @param curvesData {list} first parameter whose length should be displayed
#' @param treeData {list} second parameter whose length should be displayed
#' @param treeConsensus {list} third parameter whose length should be displayed
Experiment0.__correctnessCheck <- function(curvesData, treeData, treeConsensus) {
  print(length(curvesData));
  print(length(treeData));
  print(length(treeConsensus));
}