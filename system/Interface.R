#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

INTERFACE_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Allows to execute the implemented approaches.
#' Hint: All approaches require a gap free chronology and samples!


#########################################################
#' Executes the Consensus Approach per sample.
#'
#' @param consensusPath {string} the path to the consensus
#' @param consensusName {string} the name of the consensus file
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param bestYearsMax {numeric} tells how many best ranked years should be stored
#' @param save {logical} tells if the list of possible dates should be stored on hard disk
#' @param fileName {string} the filename without extension for the stored per sample dates-file
#'
#' @return {matrix} the matrix of possible dates
#' @examples 
#' #consensus
#' year, density, characteristic (optional - not used)
#' 1992, 0.23,    0.35
#' ...
#' 1992, 0.54,    0.35
#' 1993, 0.17,    0.52
#' ...
#' 1993, 0.67,    0.52
#' ...
#' 
#' #samples
#' part, density, characteristic (optional - not used)
#' 1,     0.17,   0.51
#' 1,     0.23,   0.51
#' ...
#' 1,     0.57,   0.51
#' 2,     0.18,   0.65
#' 2,     0.19,   0.65
#' ...
#' 2,     0.59,   0.65
#' ...
#' 
#' #score types
#' (a): Average Point-Distance
#' (b): Average Point-Distance on Slopes
#' (c): Average Point-Distance with z-scores
#' (d): Average Point-Distance on Slopes with z-scores
#' 
#' #matrix
#' sample,  rank[1],  rank[2], ..., rank[n]
#' 401,     2015,     2017,     .., 1922
#' 402,     2013,     2030,     .., 1629
#' ...
#' @export
Interface.computeDatesConsensusApproach <- function(consensusPath, consensusName, samplesPath, scoreType, 
                                                    bestYearsMax, save, fileName) {
  # load
  masterChronologyData <- Loader.readInProfiles(consensusPath, consensusName, Symbols.REAL_DATA_SEPARATOR, 
                                                columnName = Defaults.DENSITY_NAME);
  samplesData <- Loader.readInSamples(samplesPath, Symbols.REAL_DATA_SEPARATOR, Extensions.CSV);
  
  # compute scores
  scoresPerSample <- list();
  
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, samplesData,
                                                       FALSE, FALSE);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, samplesData,
                                                       FALSE, TRUE);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, samplesData,
                                                       TRUE, FALSE);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample(masterChronologyData, samplesData,
                                                       TRUE, TRUE);
  }

  # predict years
  yearsPerSample <- Analyzer.getYearsPerSample(scoresPerSample, samplesData, masterChronologyData$years, bestYearsMax, list());
  
  # store per sample predicted years
  if (save) Storer.storeTable(yearsPerSample, fileName, Extensions.CSV);
  
  return(yearsPerSample);
}


#########################################################
#' Executes the Bucket Approach per sample.
#'
#' @param bucketsPath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param innerFunc {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the given position (e.g. the minimum function min)
#' @param outerFunc {function} the function which should be applied on the 
#' per sample computed scores to get a final score for the position (e.g. the summation function sum)
#' @param bestYearsMax {numeric} tells how many best ranked years should be stored
#' @param qualityMeasures {string} the string telling you which quality measures should be active 
#' (combine multiple options: "" = none, "p" = p-values, "s" = scores, "ps" or "sp" = scores and p-values)
#' @param save {logical} tells if the list of possible dates should be stored on hard disk
#' @param fileName {string} the filename without extension for the stored per sample dates-file
#' 
#' @return {matrix} the matrix of possible dates
#' @examples 
#' #consensi for buckets
#' year, density, characteristic (optional - not used)
#' 1992, 0.23,    0.35
#' ...
#' 1992, 0.54,    0.35
#' 1993, 0.17,    0.52
#' ...
#' 1993, 0.67,    0.52
#' ...
#' 
#' (further examples: Interface.computeDatesConsensusApproach)
#' 
#' @export
Interface.computeDatesBucketApproach <- function(bucketsPath, samplesPath, 
                                                 scoreType, innerFunc, outerFunc, 
                                                 bestYearsMax, qualityMeasures = Symbols.EMPTY, 
                                                 save, fileName) {
  
  scoresData <- Interface.__getBucketsScores(bucketsPath, samplesPath, 
                                             scoreType, innerFunc, outerFunc);
  
  scoresPerSample <- scoresData$scoresPerSample;
  samplesData <- scoresData$samplesData;
  years <- scoresData$years;
  
  # interprete quality string
  qualities <- Interpreter.interpreteQualityString(qualityMeasures);
  
  # predict years
  yearsPerSample <- Analyzer.getYearsPerSample(scoresPerSample, samplesData, years, bestYearsMax, qualities);
  
  # store per sample predicted years
  if (save) Storer.storeTable(yearsPerSample, fileName, Extensions.CSV);
  
  return(yearsPerSample);
} 


#########################################################
#' Returns the buckets-scores.
#' 
#' @param bucketsPath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param innerFunc {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the position (e.g. the minimum function min)
#' @param outerFunc {function} the function which should be applied on the 
#' per sample computed scores to get a final score for the position (e.g. the summation function sum)
#'
#' @return {list} the scores per sample
#' @export
Interface.__getBucketsScores <- function(bucketsPath, samplesPath, 
                                         scoreType, innerFunc, outerFunc) {
  # load
  samplesData <- Loader.readInSamples(samplesPath, Symbols.REAL_DATA_SEPARATOR, Extensions.CSV);
  
  perYearProfiles <- new.env();
  curvesData <- Loader.readInDefaultConsensi(bucketsPath);
    
  # retrieve per year profiles
  perYearProfilesData <- CurvesMiner.collectData(curvesData);
  perYearProfiles <- perYearProfilesData$profilesPerYear;
  
  # retrieve chronology years
  startYear <- perYearProfilesData$minStartYear;
  endYear <- perYearProfilesData$maxStartYear;
  years <- startYear:endYear;

  # compute scores
  scoresPerSample <- list();
  
  # search samples in bucket chronology
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        FALSE, FALSE, innerFunc, outerFunc);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        FALSE, TRUE, innerFunc, outerFunc);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        TRUE, FALSE, innerFunc, outerFunc);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        TRUE, TRUE, innerFunc, outerFunc);
  }
  
  return(list(scoresPerSample = scoresPerSample, samplesData = samplesData, years = years));
}


#########################################################
#' Executes the Voting Approach per sample.
#' Hint: NA's are possible since it can 
#' happen that e.g. for all series-lengths
#' the same year is selected.
#'
#' @param bucketsPath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param topYearsCount {numeric} the number of years selected per column
#' @param approach {string} the string telling you which approaches should be active 
#' (multiple options - powerset approach: "p", empty string to select none)
#' @param minimumLength {numeric} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit)
#' @param bestYearsMax {numeric} tells how many best ranked years should be stored
#' @param save {logical} tells if the list of possible dates should be stored on hard disk
#' @param fileName {string} the filename without extension for the stored per sample dates-file
#'
#' @return {matrix} the matrix of possible dates
#' @export
Interface.computeDatesVotingApproach <- function(bucketsPath, samplesPath, 
                                                 scoreType, topYearsCount, approach = Symbols.EMPTY, 
                                                 minimumLength, bestYearsMax, save, fileName) {
  
  scoresData <- Interface.__getBucketsScores(bucketsPath, samplesPath, 
                                             scoreType, min, identity);
  
  scoresPerSample <- scoresData$scoresPerSample;
  samplesData <- scoresData$samplesData;
  years <- scoresData$years;
  
  # interprete quality string
  qualities <- Interpreter.interpreteApproachString(approach);
  
  # create scores table
  scoresTables <- Analyzer.getScoresTables(scoresPerSample, years);
 
  # predict years
  topYearsPerSample <- Analyzer.getTopYearsPerSample(scoresTables, topYearsCount, 
                                                     qualities$powersetApproach, qualities$doubleWeighting, 
                                                     qualities$logarithmicWeighting, FALSE, minimumLength);

  yearsPerSample <- Analyzer.getYearsPerSample2(topYearsPerSample, samplesData, bestYearsMax);
  
  # store per sample predicted years
  if (save) Storer.storeTable(yearsPerSample, fileName, Extensions.CSV); 
  
  return(yearsPerSample);
}


#########################################################
#' Executes first a fast Point-Based Approach (correlation coefficient / t-value based approach)
#' and afterwards the Bucket Approach on an amount of years found by the Points-Based Approach.
#' Hint: topYearsCount has to be set at least to 30 to get a proper distribution for p-values!
#'
#' @param bucketsPath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreTypeCharacteristic {string} the score-type which should be used for computation during the characteristic approach 
#' ("p" = Pearson's Rho, "t" = Kendall's Tau, "r" =  Spearman's Rho, "v" = t-value)
#' @param scoreTypeBucket {string} the score-type which should be used for computation during the buckets-approach
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param topYearsCount {numeric} the number of top years which are stored by the characteristic approach (at least 2)
#' @param bestYearsMax {numeric} tells how many best ranked years should be stored
#' @param qualityMeasures {string} the string telling you which quality measures should be active 
#' (combine multiple options: "" = none, "p" = p-values, "s" = scores, "ps" or "sp" = scores and p-values)
#' @param save {logical} tells if the list of possible dates should be stored on harddisk
#' @param fileName {string} the filename without extension for the stored per sample dates-file
#'
#' @return {matrix} the matrix of possible dates
#' @export
Interface.computeDatesTwoStepApproach <- function(bucketsPath, samplesPath, 
                                                  scoreTypeCharacteristic, scoreTypeBucket,
                                                  topYearsCount, bestYearsMax, qualityMeasures, 
                                                  save, fileName) {
  
  if (topYearsCount == 1) stop(Exceptions.TOP_YEARS);
  
  # read in data
  samplesData <- Loader.readInSamples(samplesPath, Symbols.REAL_DATA_SEPARATOR, Extensions.CSV);
  curvesData <- Loader.readInDefaultConsensi(bucketsPath);
  
  # retrieve data
  # profiles
  perYearProfilesData <- CurvesMiner.collectData(curvesData);
  perYearProfiles <- perYearProfilesData$profilesPerYear;
  
  startYear <- perYearProfilesData$minStartYear;
  endYear <- perYearProfilesData$maxStartYear;
  years <- startYear:endYear;
  
  # compute characteristic scores
  characteristicScoresPerSample <- Interface.__getRingCharacteristicScores(curvesData, samplesData, 
                                                                           scoreTypeCharacteristic, 
                                                                           startYear, endYear);

  # create scores table
  scoresTables <- Analyzer.getScoresTables(characteristicScoresPerSample, years);
  
  # store top rated positions or years with the highest correlation coefficients
  topYearsSample <- Analyzer.getTopYearsPerSample(scoresTables, topYearsCount, FALSE, 
                                                  FALSE, FALSE, TRUE);
  
  yearsIndices <- Analyzer.getIndicesPerYear(years);
  
  # compute scores only for given years in master-chronology
  scoresPerSample <- list();
  
  # search samples in bucket chronology
  if (scoreTypeBucket == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        FALSE, FALSE, min, sum, 
                                                        topYearsSample, yearsIndices);
  } else if (scoreTypeBucket == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        FALSE, TRUE, min, sum,
                                                        topYearsSample, yearsIndices);
  } else if (scoreTypeBucket == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        TRUE, FALSE, min, sum,
                                                        topYearsSample, yearsIndices);
  } else if (scoreTypeBucket == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample2(perYearProfiles, samplesData,
                                                        TRUE, TRUE, min, sum,
                                                        topYearsSample, yearsIndices);
  }
  
  # predict years
  qualities <- Interpreter.interpreteQualityString(qualityMeasures);

  yearsPerSample <- Analyzer.getYearsPerSample3(scoresPerSample, samplesData, 
                                                topYearsSample, bestYearsMax, qualities);

  # store per sample predicted years
  if (save) Storer.storeTable(yearsPerSample, fileName, Extensions.CSV); 
  
  return(yearsPerSample);
}


#########################################################
#' Returns the scores for characteristics and maximum densities.
#'
#' @param curvesData {list} data of multiple curves
#' @param samplesData {list} data of samples
#' @param scoreType {string} the score-type which should be used for computation 
#' ("p" = Pearson's Rho, "t" = Kendall's Tau, "r" =  Spearman's Rho, "v" = t-value)
#' @param startYear {numeric} the year in which the characteristic chronology should start
#' @param endYear {numeric} the year in which the characteristic chronology should end
#'
#' @return {list(scoresPerSample, samplesData, perYearProfiles, years)} the scores per sample
#' @export
Interface.__getRingCharacteristicScores <- function(curvesData, samplesData, scoreType, startYear, endYear) {

  # characteristics - due to limited amount of time characteristics are renamed as "widths" (to avoid earlier code and breaking)
  perYearCharacteristicsData <- CurvesMiner.collectWidthsData(curvesData);  
  perYearCharacteristics <- perYearCharacteristicsData$widthsPerYear;
  
  masterChronology <- CurvesMiner.computeCharacteristicsChronology(perYearCharacteristics, startYear, endYear);
  
  scoresPerSample <- list();
  
  if (scoreType == Defaults.TASK_PEARSON_CHAR) {
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(masterChronology$widths, samplesData, 
                                                             Defaults.CORRELATION_PEARSON, FALSE);
  } else if (scoreType == Defaults.TASK_KENDALL_CHAR) {
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(masterChronology$widths, samplesData, 
                                                             Defaults.CORRELATION_KENDALL, FALSE);
  } else if (scoreType == Defaults.TASK_SPEARMAN_CHAR) {
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(masterChronology$widths, samplesData, 
                                                             Defaults.CORRELATION_SPEARMAN, FALSE);
  } else if (scoreType == Defaults.TASK_T_VALUE_CHAR) {
    scoresPerSample <- Analyzer.computeCoefficientsPerSample(masterChronology$widths, samplesData, 
                                                             Defaults.CORRELATION_T_VALUE, FALSE);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Executes the Bucket Approach per sample.
#'
#' @param perTreePath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#' @param bestYearsMax {numeric} tells how many best ranked years should be stored
#' @param save {logical} tells if the list of possible dates should be stored on hard disk
#' @param fileName {string} the filename without extension for the stored per sample dates-file
#'
#' @return {matrix} the matrix of possible dates 
#' @export
Interface.computeDatesPerTreeApproach <- function(perTreePath, samplesPath, 
                                                  scoreType, bestYearsMax, save, fileName) {
  scoresData <- Interface.__getPerTreeScores(perTreePath, samplesPath, scoreType);
  
  scoresPerSample <- scoresData$scoresPerSample;
  samplesData <- scoresData$samplesData;
  years <- scoresData$years;

  # predict years
  yearsPerSample <- Analyzer.getYearsPerSample(scoresPerSample, samplesData, years, 
                                               bestYearsMax, list());
  
  # store per sample predicted years
  if (save) Storer.storeTable(yearsPerSample, fileName, Extensions.CSV);
  
  return(yearsPerSample);
}


#########################################################
#' Returns the per-tree computed scores.
#'
#' @param perTreePath {string} the path to the per tree consensi which should be used
#' @param samplesPath {string} the path to the samples which should be dated
#' @param scoreType {string} the score-type which should be used for computation 
#' ("a" = y-based, "b" = slope-based, "c" = z-scores y-based, "d" = z-scores slope-based)
#'
#' @return {list} the scores per sample
#' @export
Interface.__getPerTreeScores <- function(perTreePath, samplesPath, scoreType) {
  # load
  samplesData <- Loader.readInSamples(samplesPath, Symbols.REAL_DATA_SEPARATOR, Extensions.CSV);
  curvesData <- Loader.readInDefaultConsensi(perTreePath);
  
  # retrieve chronology years
  yearsData <- Analyzer.getStartAndEndYears(curvesData);
  startYears <- yearsData$startYears;
  endYears <- yearsData$endYears;
  
  collectedData <- CurvesMiner.collectData(curvesData);
  minStartYear <- as.numeric(collectedData$minStartYear);
  maxEndYear <- as.numeric(collectedData$maxStartYear);
  
  # compute scores
  scoresPerSample <- list();
  
  # search samples in bucket chronology
  if (scoreType == Defaults.TASK_A) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(curvesData, samplesData,
                                                        FALSE, FALSE, startYears, endYears, 
                                                        minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_B) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(curvesData, samplesData,
                                                        FALSE, TRUE, startYears, endYears, 
                                                        minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_C) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(curvesData, samplesData,
                                                        TRUE, FALSE, startYears, endYears, 
                                                        minStartYear:maxEndYear);
  } else if (scoreType == Defaults.TASK_D) {
    scoresPerSample <- Analyzer.computeScoresPerSample4(curvesData, samplesData,
                                                        TRUE, TRUE, startYears, endYears, 
                                                        minStartYear:maxEndYear);
  }
  
  return(list(scoresPerSample = scoresPerSample, samplesData = samplesData, years = minStartYear:maxEndYear));
}