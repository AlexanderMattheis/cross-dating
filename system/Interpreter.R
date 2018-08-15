#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

INTERPRETER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Allows to interpret commands from the Interface.


#########################################################
#' Interprets the quality string and returns options, 
#' conretly p-values activation state and the scores activation state dependent on the string.
#'
#' @param qualityMeasures {string} the string which has to be interpreted
#'
#' @return {list(pValues, scores)} the booleans telling you, if the given
Interpreter.interpreteQualityString <- function(qualityMeasures) {
  # divide string
  options <- unlist(stringr::str_split(qualityMeasures, Symbols.EMPTY));
  
  pValues <- FALSE;
  scores <- FALSE;
  
  if (length(options) > 0) {
    for (i in 1:length(options)) {
      option <- options[[i]];
      
      if (option == Defaults.P_VALUES_COMMAND) {
        pValues <- TRUE;
      } else if (option == Defaults.SCORES_COMMAND) {
        scores <- TRUE;
      }
    }
  } 

  return(list(pValues = pValues, scores = scores));
}


#########################################################
#' Interprets the approach string and activates the corresponding options.
#'
#' @param approaches {string} the string which has to be interpreted
#'
#' @return {list(doubleWeighting, logarithmicWeighting, powersetApproach)} 
#' the states telling you if given approaches are activated or deactivated
Interpreter.interpreteApproachString <- function(approaches) {
  # divide string
  options <- unlist(stringr::str_split(approaches, Symbols.EMPTY));
  
  doubleWeighting <- FALSE;
  logarithmicWeighting <- FALSE;
  powersetApproach <- FALSE;
  
  if (length(options) > 0) {
    for (i in 1:length(options)) {
      option <- options[[i]];
      
      if (option == Defaults.DOUBLE_WEIGHTING_COMMAND) {
        doubleWeighting <- TRUE;
      } else if (option == Defaults.LOG_WEIGHTING_COMMAND) {
        logarithmicWeighting <- TRUE;
      } else if (option == Defaults.POWER_SET_APPROACH_COMMAND) {
        powersetApproach <- TRUE;
      }
    }
  } 
  
  return(list(doubleWeighting = doubleWeighting, 
              logarithmicWeighting = logarithmicWeighting,
              powersetApproach = powersetApproach));
}