#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

ANALYZER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Contains functions to analyze and edit data.

#########################################################
#' Every subpattern starting at patterns first position/profile is aligned 
#' against every consensus-part to get a distance.
#' These distances are then returned.
#'
#' @param consensusProfiles {list} the list of profiles of the consensus
#' @param patternProfiles {list} the list of profiles of the pattern
#' @param distanceFunction {numerical} the distance function which should be used to measure distances between profiles
#' 
#' @return {list} the scores per pattern length over all possible positions
#' @export
Analyzer.computeMicaScores <- function(consensusProfiles, patternProfiles, distanceFunction) {
  scoreListsPerPatternLength <- list();
  
  n <- length(consensusProfiles);
  m <- length(patternProfiles);
  
  # for pattern/profile series of ascending length
  for (i in 1:m) {
    scores <- Analyzer.__computePositionWiseScores(n, i, consensusProfiles, patternProfiles, distanceFunction);
    scoreListsPerPatternLength <- rlist::list.append(scoreListsPerPatternLength, scores);
  }
  
  return(scoreListsPerPatternLength);
}


#########################################################
#' Computes the score for each position of the pattern in the consensus. 
#' Hint: With this function, 
#' it is possible to determine 
#' how many profiles are needed to date a pattern.
#'
#' @param n {numerical} number of profiles in consensus curve
#' @param m {numerical} number of profiles in pattern curve
#' @param consensusProfiles {list} the list of profiles of the consensus
#' @param patternProfiles {list} the list of profiles of the pattern
#' @param distanceFunction {numerical} the distance function which should be used to measure distances between profiles
#'
#' @return {list} the scores for each position of the pattern in the consensus
Analyzer.__computePositionWiseScores <- function(n, m, consensusProfiles, patternProfiles, distanceFunction) {
  scores <- list()
  
  # iterate over each curve position with the subpattern
  for (j in 1:(n-m+1)) {  # "+1", because i is "1" at the beginning and else we wouldn't align against every position 
    subconsensusProfiles <- CurvesMiner.getSubset(consensusProfiles, j, (j-1) + m);  # "-1", because j is "1" at the beginning
    subpatternProfiles <- CurvesMiner.getSubset(patternProfiles, 1, m);
    
    accumulatedScore <- 0;
    
    # iterate over each current position of the pattern in consensus curve 
    # and align the two profiles: one from consensus, one from subpattern and sum up the scores
    for (k in 1:m) {
      consensusProfile <- subconsensusProfiles[[k]];
      subpatternProfile <- subpatternProfiles[[k]];
      
      # compute score
      profilesToAlign <- CurvesMiner.__equalizeAllPointNumbers(list(consensusProfile, subpatternProfile));  # encoding
      alignmentScore <- Analyzer.__computeScore(profilesToAlign, distanceFunction);
      
      accumulatedScore <- accumulatedScore + alignmentScore;
    }
    
    scores <- rlist::list.append(scores, accumulatedScore);  
  }
  
  return(scores);
}


#########################################################
#' Aligns the profiles in the given dataframe and returns their aligned distance.
#' Hint: First profile is the reference, which is not warped.
#' Hint 2: Sadly to use pairwise distance is wrong! So, this function is only for testing.
#'
#' @param profilesToAlign {dataframe} the two profiles, one from pattern and one from curve
#' @param distanceFunction {numerical} the distance function
#'
#' @return {numerical} the score between aligned curve and pattern profile
Analyzer.__computeScore <- function(profilesToAlign, distanceFunction) {
  output <- alignCurves(y = profilesToAlign, 
                        distFunc = distanceFunction,
                        distSample = Parameters.DIST_SAMPLE,
                        distWarpScaling = Parameters.DIST_WARP_SCALING,
                        maxWarpingFactor = Parameters.MAX_WARPING_FACTOR,
                        maxRelXShift = Parameters.MAX_REL_X_SHIFT,
                        minRelIntervalLength = Parameters.MIN_REL_INTERVAL_LENGTH,
                        minRelMinMaxDist = Parameters.MIN_REL_MIN_MAX_DIST,
                        minRelSlopeHeight = Parameters.MIN_REL_SLOPE_HEIGHT,
                        reference = 1,  # to use the first profile as reference
                        outSlope = Parameters.OUT_SLOPE);
  
  alignmentScore <- output$pairDist$warped[1, 2];  # matrix of pairwise scores at position [1, 2] (dist. between profile 1 and 2)
  return(alignmentScore);
}


#########################################################
#' Interpolates profiles to the same length
#' and then computes the mean distance between sample
#' and current consensus cut. This is done for each
#' position of the consensus.
#' Hint: Tasks (D) a)-d).
#'
#' @param consensusProfiles {list} the list of profiles of the consensus (profiles with its years)
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#' 
#' @return {list} the scores per sample
Analyzer.computeScoresPerSample <- function(consensusProfiles, samples,
                                            normalize, slopeBased) {
  scoresPerSample <- list();
  
  # prepare data
  data <- Analyzer.__dataPreparation(consensusProfiles, samples, normalize, slopeBased);
  
  print(Titles.CURRENT_PROCESSED_SAMPLE);
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    patternProfiles <- data$samples[[i]]$profiles;
    scores <- Analyzer.__computeScores(data$consensusProfiles, patternProfiles, FALSE, slopeBased);
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Prepares the data i.e. replaces with the slopes or normalized data.
#'
#' @param consensusProfiles {list} the list of profiles of the consensus 
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#'
#' @return {list(consensusProfiles, samples)}
Analyzer.__dataPreparation <- function(consensusProfiles, samples, normalize, slopeBased) {
  consensusProfiles <- Math.computeOperation(list(consensusProfiles), normalize, slopeBased);
  samples <- Math.computeOperation(samples, normalize, slopeBased);
  
  return(list(consensusProfiles = consensusProfiles[[1]][[1]], samples = samples));
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the consensus, which consists out of profiles.
#'
#' @param consensusProfiles {list} the list of profiles of the consensus
#' @param patternProfiles {list} the list of profiles of the pattern
#' @param mica {logical} tells if mica should be used or not
#' @param normalized {logical} tells if mica is set to true, 
#' so if the values are normalized and the x-warping should be applied on real coordinates
#' @param slopeBased {logical} tells if mica is set to true, 
#' if the scores should be computed on the slopes
#' 
#' @return {list} the scores
Analyzer.__computeScores <- function(consensusProfiles, patternProfiles, mica=FALSE, slopeBased=FALSE) {
  scores <- list();
  
  n <- length(consensusProfiles);
  m <- length(patternProfiles);
  
  # iterate over each curve position with the subpattern
  for (i in 1:(n-m+1)) {  # "+1", because i is "1" at the beginning and else we wouldn't align against every position 
    # grab a consensus-cut of pattern-size
    consensusCut <- CurvesMiner.getSubset(consensusProfiles, i, (i-1) + m);  # "-1", because i is "1" at the beginning
    
    accumulatedScore <- 0;
    
    # iterate over each position of the pattern and the consensus-cut
    for (k in 1:m) {
      consensusProfile <- consensusCut[[k]];
      subpatternProfile <- patternProfiles[[k]];
      
      if (!mica) {
        alignmentScore <- Analyzer.__getScore(subpatternProfile, consensusProfile);
      } else {
        alignmentScore <- Analyzer.__getMicaScore(subpatternProfile, consensusProfile, slopeBased);
      }
      
      accumulatedScore <- accumulatedScore + alignmentScore;
    }
    
    #print(accumulatedScore);
    scores <- rlist::list.append(scores, accumulatedScore);  
  }
  
  return(scores);
}


#########################################################
#' Computes a score with the given parameters 
#' for the two profiles.
#'
#' @param subpatternProfile {list} the y-values of the profile
#' @param consensusProfile {list} the y-values of the profile
#'
#' @return {numerical} the score with the given parameters and properties
Analyzer.__getScore <- function(subpatternProfile, consensusProfile) {
  
  # interpolation of x-values is necessary i.e. both curves should have the same number
  # Hint: does exactly the same like MICA interpolation function, but significantly faster
  consensusProfile <- Math.interpolatePoints(consensusProfile, Parameters.RANKING_DIST_SAMPLE);
  subpatternProfile <- Math.interpolatePoints(subpatternProfile, Parameters.RANKING_DIST_SAMPLE);
  
  # compute score
  alignmentScore <- Math.getMeanDistance(consensusProfile, subpatternProfile);
  
  return(alignmentScore);
}


#########################################################
#' Does x-Warping with MICA and interpolation to the same length
#' and then computes the mean distance between sample
#' and current consensus cut. This is done for each
#' position of the consensus.
#' Hint: Tasks (D) e)-h).
#'
#' @param consensusProfiles {list} the list of profiles of the consensus (profiles with its years)
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#' 
#' @return {list} the scores per sample
#' @export
Analyzer.computeMicaScoresPerSample <- function(consensusProfiles, samples,
                                                normalize, slopeBased) {
  scoresPerSample <- list();
  
  # prepare data
  data <- Analyzer.__dataPreparation(consensusProfiles, samples, normalize, FALSE);
  
  print(Titles.CURRENT_PROCESSED_SAMPLE);
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    patternProfiles <- data$samples[[i]]$profiles;
    scores <- Analyzer.__computeScores(data$consensusProfiles, patternProfiles, TRUE, slopeBased);
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Computes a score with the given parameters 
#' for the two profiles.
#'
#' @param subpatternProfile {list} the y-values of the profile
#' @param consensusProfile {list} the y-values of the profile
#' @param slopeBased {logical} tells if the scores should be computed on the slopes
#'
#' @return {numerical} the score with the given parameters and properties
Analyzer.__getMicaScore <- function(subpatternProfile, consensusProfile, slopeBased) {
  
  # get interpolated warped coordinate values
  profilesToAlign <- CurvesMiner.__equalizeAllPointNumbers(list(consensusProfile, subpatternProfile));
  yCurves <- Analyzer.__getAfterAlignmentYCoordinates(profilesToAlign, slopeBased);
  
  # compute score
  alignmentScore <- Math.getMeanDistance(yCurves$curve1, yCurves$curve2);
  
  return(alignmentScore);
}


#########################################################
#' Computes an alignment with MICA and gets the y-coordinates,
#' which are used to compute the scores.
#'
#' @param profilesToAlign {dataframe} the two profiles, one from pattern and one from curve
#' @param slopes {logical} if selected slopes, instead of y-values returned
#'
#' @return {list(curve1, curve2)} coordinates of all curves
Analyzer.__getAfterAlignmentYCoordinates <- function(profilesToAlign, slopes) {
  interpolatedXy <- data.matrix;
  
  if (slopes) {
    output <- alignCurves(y = profilesToAlign, 
                          distFunc = Defaults.CURVE_SLOPE_MEAN_ABSOLUTE_DISTANCE,
                          distSample = Parameters.RANKING_DIST_SAMPLE,
                          distWarpScaling = Parameters.RANKING_DIST_WARP_SCALING,
                          maxWarpingFactor = Parameters.RANKING_MAX_WARPING_FACTOR,
                          maxRelXShift = Parameters.RANKING_MAX_REL_X_SHIFT,
                          minRelIntervalLength = Parameters.RANKING_MIN_REL_INTERVAL_LENGTH,
                          minRelMinMaxDist = Parameters.RANKING_MIN_REL_MIN_MAX_DIST,
                          minRelSlopeHeight = Parameters.RANKING_MIN_REL_SLOPE_HEIGHT,
                          reference = Parameters.RANKING_REFERENCE,  # to use the first profile as reference
                          outSlope = FALSE);
    
    interpolatedXy <- getSlopeOfCurves(output$xWarped, profilesToAlign);
    interpolatedXy <- interpolateCurves(interpolatedXy$x, interpolatedXy$y, Parameters.RANKING_DIST_SAMPLE); 
  } else {  # y-values
    output <- alignCurves(y = profilesToAlign, 
                          distFunc = Defaults.CURVE_SLOPE_MEAN_ABSOLUTE_DISTANCE,
                          distSample = Parameters.RANKING_DIST_SAMPLE,
                          distWarpScaling = Parameters.RANKING_DIST_WARP_SCALING,
                          maxWarpingFactor = Parameters.RANKING_MAX_WARPING_FACTOR,
                          maxRelXShift = Parameters.RANKING_MAX_REL_X_SHIFT,
                          minRelIntervalLength = Parameters.RANKING_MIN_REL_INTERVAL_LENGTH,
                          minRelMinMaxDist = Parameters.RANKING_MIN_REL_MIN_MAX_DIST,
                          minRelSlopeHeight = Parameters.RANKING_MIN_REL_SLOPE_HEIGHT,
                          reference = Parameters.RANKING_REFERENCE,  # to use the first profile as reference
                          outSlope = FALSE);
    
    # a factor is multiplied to the length after interpolation, to avoid an information loss
    interpolatedXy <- interpolateCurves(output$xWarped, profilesToAlign, Parameters.RANKING_DIST_SAMPLE);  # interpolate y-values
  }
  
  return(list(curve1 = as.list(interpolatedXy$y[,1]), curve2 = as.list(interpolatedXy$y[,2])));
}


#########################################################
#' Extracts the score of each sample at its correct position in the consensus.
#'
#' @param scoresPerSample {list} the scores per sample
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param years {list} the years of the consensus
#'
#' @return {list} score per sample
Analyzer.getScoresForCorrectSamplePositions <- function(scoresPerSample, samples, years) {
  scores <- list();
  
  # get for each sample the correct score
  for (i in 1:length(samples)) {
    scoresOfSample <- scoresPerSample[[i]];
    sample <- samples[[i]];
    correctScore <- Analyzer.__getCorrectPatternPositionScore(scoresOfSample, 
                                                              sample$years[[1]], years);
    scores <- rlist::list.append(scores, correctScore);
  }
  
  # check up (there should not be any -1)
  print(scores);
  
  return(scores);
}


#########################################################
#' Returns the score for the pattern at its correct position.
#'
#' @param scoresOfSample {list} the scores for each year between the chronology and a sample
#' @param sampleYear {numerical} the year from which the score should be extracted
#' @param years {list} the years of the consensus (hint: gaps between years possible)
#'
#' @return {list} the scores for the given position
Analyzer.__getCorrectPatternPositionScore <- function(scoresOfSample, sampleYear, years) {
  return(Analyzer.__getCorrectPatternPositionData(scoresOfSample, sampleYear, years));
}


#########################################################
#' Returns the data for the pattern at its correct position.
#'
#' @param dataOfSample {list} the data for each year between the chronology and a sample
#' @param sampleYear {numerical} the year from which the score should be extracted
#' @param years {list} the years of the consensus (hint: gaps between years possible)
#'
#' @return {list} the data for the given position
Analyzer.__getCorrectPatternPositionData <- function(dataOfSample, sampleYear, years) {
  # iterate over each year until the sampleYear is found
  for (i in 1:length(years)) {
    year <- years[[i]];
    
    if (year == sampleYear) {  # if found
      score <- dataOfSample[[i]];  # score from this year (and consecutive, which the pattern overlaps)
      return(score);
    }
  }
  
  return(-1);
}


#########################################################
#' Returns the ranks of the given scores inbetween all scores.
#' 
#' @param scoresForCorrectPositions {list} the scores for the correct position of the pattern
#' @param scoresPerSample {list} the scores per sample
#' @param decreasing {logical} tells if the highest value should have lowest rank
#' 
#' @return {list} the ranks for the given position
Analyzer.getRanks <- function(scoresForCorrectPositions, scoresPerSample, decreasing = FALSE) {
  ranks <- list();
  
  # iterate over all correct position scores
  for (i in 1:length(scoresForCorrectPositions)) {
    # find out position of correct score in all scores of that sample
    allScores <- unlist(scoresPerSample[[i]]);
    sortedScores <- sort(allScores, decreasing = decreasing);
    correctScore <- scoresForCorrectPositions[[i]];
    ranksWithThatScore <- which(sortedScores == correctScore);  # position of correct score
    rankOfCorrectScore <- ranksWithThatScore[[1]] + length(ranksWithThatScore) - 1;  # last one
    ranks <- rlist::list.append(ranks, rankOfCorrectScore);
  }
  
  return(ranks);
}


#########################################################
#' Interpolates profiles to the same length
#' and then computes the mean distance between a sample
#' and a bucket. This is done for each
#' position of the chronology.
#' Hint: Tasks (D) a)-d).
#'
#' @param bucketChronology {list} the list of buckets
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#' @param func {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the position
#' @param func2 {function} the function which should be applied on the 
#' per sample computed scores to get a score for the position
#' @param topYearsPerSample {list} the years at which it should be looked at per sample
#' @param yearsIndices {environment(years, indices)} the indices per year
#' 
#' @return {list} the scores per sample
Analyzer.computeScoresPerSample2 <- function(bucketChronology, samples,
                                             normalize, slopeBased, func, func2, 
                                             topYearsPerSample = list(),
                                             yearsIndices = env()) {
  scoresPerSample <- list();
  
  # prepare data
  data <- Analyzer.__dataPreparation2(bucketChronology, samples, normalize, slopeBased);
  
  print(Titles.CURRENT_PROCESSED_SAMPLE);
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    patternProfiles <- data$samples[[i]]$profiles;
    
    if (length(topYearsPerSample) > 0) {
      topYears <- unlist(topYearsPerSample[[i]]);
      scores <- Analyzer.__computeBucketScores(data$bucketChronology, patternProfiles, func, func2, 
                                               topYears, yearsIndices);
    } else {
      scores <- Analyzer.__computeBucketScores(data$bucketChronology, patternProfiles, func, func2);
    }
    
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Prepares the data i.e. replaces with the slopes or normalized data.
#'
#' @param bucketChronology {envrionment} a chronology of buckets
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute a score on slopes
#'
#' @return {list(bucketChronology, samples)} the prepared data
Analyzer.__dataPreparation2 <- function(bucketChronology, samples, normalize, slopeBased) {
  editedBucketChronology <- list();
  
  # iterate over all years and compute for every profile in a year the desired operations
  for (i in ls(bucketChronology)) {
    bucket <- bucketChronology[[i]];
    print(i);
    
    if (!is.null(bucket)) {  # if bucket exist
      editedBucket <- list();
      
      # iterate over all profiles in bucket
      for (j in 1:length(bucket)) {
        data <- bucket[[j]];
        
        if (normalize) data <- Math.computePureValues(data);
        if (slopeBased) data <- Math.computeDerivate(data);
        
        editedBucket <- rlist::list.append(editedBucket, data);
      }
      
      editedBucketChronology <- rlist::list.append(editedBucketChronology, editedBucket);
    }
  }
  
  samples <- Math.computeOperation(samples, normalize, slopeBased);
  
  return(list(bucketChronology = editedBucketChronology, samples = samples));
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the bucket-chronology, which consists out of buckets of profiles.
#'
#' @param bucketChronology {list} the list of buckets
#' @param patternProfiles {list} the list of profiles of the pattern
#' @param func {function} the function which should be applied on the 
#' per bucket computed scores to get a score for the position
#' @param func2 {function} the function which should be applied on the 
#' per sample computed scores to get a score for the position
#' @param topYears {list} the years at which it should be looked at
#' @param yearsIndices {environment(years, indices)} the indices per year
#' 
#' @return {list} the scores
Analyzer.__computeBucketScores <- function(bucketChronology, patternProfiles, func, func2, 
                                           topYears = list(), yearsIndices = list()) {
  scores <- list();
  
  n <- length(bucketChronology);
  m <- length(patternProfiles);
  
  indicesToLookAt <- 1:(n-m+1);
  
  if (length(topYears) > 0) {  # extension for ring-width based heuristic
    # find the right indices
    indicesToLookAt <- list();
    
    for (year in topYears) {
      index <- yearsIndices[[toString(year)]];
      indicesToLookAt <- rlist::list.append(indicesToLookAt, index);
    }
    
    indicesToLookAt <- unlist(indicesToLookAt);
  }
  
  # iterate over each curve position with the subpattern
  for (i in indicesToLookAt) {  # "+1", because i is "1" at the beginning and else we wouldn't align against every position
    # grab a cut of pattern-size
    bucketCut <- CurvesMiner.getSubset(bucketChronology, i, (i-1) + m);  # "-1", because i is "1" at the beginning
    
    accumulatedScore <- 0;
    scoresPerSampleProfile <- list();
    
    # iterate over each position of the pattern and the bucket-cut
    for (k in 1:m) {
      bucket <- bucketCut[[k]];
      subpatternProfile <- patternProfiles[[k]];
      
      bucketScores <- list();
      # print(length(bucket));
      
      # iterate over all profiles in the bucket
      for (l in 1:length(bucket)) {
        bucketProfile <- bucket[[l]];
        
        if (identical(subpatternProfile, bucketProfile)) { 
          print(Strings.BUG);  # this should be impossible
          View(subpatternProfile);
          View(bucketChronology);
          stop("DONE");
        }
        
        alignmentScore <- Analyzer.__getScore(subpatternProfile, bucketProfile);
        
        bucketScores <- rlist::list.append(bucketScores, alignmentScore);
      }
      
      scoresPerSampleProfile <- rlist::list.append(scoresPerSampleProfile, func(unlist(bucketScores)));
    }
    
    # add together
    aggregatedScores <- unlist(scoresPerSampleProfile);
  
    scores <- rlist::list.append(scores, func2(aggregatedScores)); 
  }
  
  return(scores);
} 


#########################################################
#' Returns the pairwise distances of the profile to the profiles in the list.
#'
#' @param profilesList {list} the list from which profiles are picked out
#' @param profile {list} the profile for which the distances are computed
#' @param normalize {logical} tells if the profiles should be normalized before computing any distances
#'
#' @return {list(pairwiseDistance)} the list of pairwise distances for the list
Analyzer.getPairwiseDistances <- function(profilesList, profile, normalize) {
  pairwiseDistances <- list();
  
  if (length(profilesList) > 0) {
    # iterate over all profiles in the list
    for (i in 1:length(profilesList)) {
      profileFromList <- profilesList[[i]];
      
      if (normalize) {
        profileFromList <- Math.computePureValues(profileFromList);
        profile <- Math.computePureValues(profile);
      }
      
      distance <- Math.getMeanDistance(profileFromList, profile);
      pairwiseDistances <- rlist::list.append(pairwiseDistances, distance);
    }
  }
  
  return(pairwiseDistances);
}


#########################################################
#' Returns the top years per sample.
#'
#' @param scoresTable {list} the score-dataframes of multiple samples with its years in the first column
#' @param topYears {numerical} the number of years to return
#' @param powersetApproach {logical} tells if the powerset approach have to be used 
#' i.e. computing a score for every subset (hint: top-years count does not used anymore, if activated)
#' @param doubleWeighting {logical} tells if the years with ranks > 1 should be weighted down
#' @param logarithmicWeighting {logical} tells if the weights should be logarimized and rounded
#' @param descending {logical} tells how to sort the scores
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
#' 
#' @return {list} the lists containing best years per column
#' @export
Analyzer.getTopYearsPerSample <- function(scoresTables, topYears, powersetApproach, 
                                          doubleWeighting = FALSE, logarithmicWeighting = FALSE, 
                                          descending = FALSE, minimumLength = -1) {
  bestYearsPerSample <- list();

  print(Titles.CURRENTLY_PROCESSED);
  
  # iterate over each table
  for (i in 1:length(scoresTables)) {
    print(i);
    scoreTable <- scoresTables[[i]];
    bestYearsPerColumn <- Analyzer.getTopYearsPerColumn(scoreTable, topYears, powersetApproach, 
                                                        doubleWeighting, logarithmicWeighting,
                                                        descending, minimumLength);
    bestYearsPerSample <- rlist::list.append(bestYearsPerSample, bestYearsPerColumn);
  }
  
  return(bestYearsPerSample);
}


#########################################################
#' Returns the top years per column.
#'
#' @param scoresTable {dataframe} the scores of a single sample with the years in the first column
#' @param topYears {numerical} the number of years to return
#' @param powersetApproach {logical} tells if the powerset approach have to be used 
#' i.e. computing a score for every subset (hint: top-years count does not used anymore, if activated)
#' @param doubleWeighting {logical} tells if the years with ranks > 1 should be weighted down
#' @param logarithmicWeighting {logical} tells if the weights should be logarimized and rounded
#' @param descending {logical} tells how to sort the scores
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
#' 
#' @return {list} the best years for each column
#' @export
Analyzer.getTopYearsPerColumn <- function(scoresTable, topYears, powersetApproach, 
                                          doubleWeighting = FALSE, logarithmicWeighting = FALSE,
                                          descending = FALSE, minimumLength = -1) {
  # test for powerset approach
  if (powersetApproach) { 
    powersetTableData <- Analyzer.__getPowerSetTable(scoresTable, minimumLength); 
    
    scoresTable <- powersetTableData$scoresTable;
    combinations <- powersetTableData$combinations;
  }

  bestYearsPerColumn <- list();
  years <- scoresTable$year;
  
  if (nrow(scoresTable) >= topYears) {  # to avoid choosing to many
    # iterate over each column and sort the scores within together with the years
    for (i in 2:ncol(scoresTable)) {
      scores <- scoresTable[i];
      scoresOrder <- order(scores, decreasing = descending);
      
      sortedYears <- years[scoresOrder];
      bestYears <- sortedYears[1:topYears];
      
      if (powersetApproach) {
        combination <- combinations[[i-1]];
        
        if (logarithmicWeighting) {
          weighting <- round(log(length(combination)));  # natural logarithm
        } else {
          weighting <- length(combination);
        }

        bestYearsPerColumn <- Analyzer.__getWeightedYearsPerColumn(bestYearsPerColumn, weighting, 
                                                                   bestYears, doubleWeighting);
      } else {
        bestYearsPerColumn <- rlist::list.append(bestYearsPerColumn, bestYears);
      }
    }
  }
  
  return(bestYearsPerColumn);
}


#########################################################
#' Computes out of the score table the powerset table.
#'
#' @param scoresTable {dataframe} the score table for which the powerset table should be computed
#' @param minimumLength {numerical} tells which minimum sample lengths should be considered
#' in the powerset table (if it is -1, then there is no limit) -> decreases the runtime
#' 
#' @return {list(scoresTable, combinations)} extended score table with the names of the combinations
#'
#' @examples 
#'       {a, b, c} -> {a}, {b}, {c}, {a + b}, {a + c}, {b + c}, {a + b + c}
#' 1922   1  2  3  ->  1    2    3      3        4        5         6  
#' 1923   4  5  6  ->  4    5    6      9      ....
#' ...
Analyzer.__getPowerSetTable <- function(scoresTable, minimumLength = -1) {
  # to remove unnecessary combinations
  start <- 1;
  
  if (minimumLength > 1) {
    start <- minimumLength;
  }
  
  # create all column combinations
  # idea taken from https://stackoverflow.com/questions/49570793/r-list-all-combinations-with-combn-multiple-m-values
  columns <- ncol(scoresTable) - 1;
  combinations <- vector("list", columns);
  
  if (start == 1) {
    firstCombi <- data.frame(X1 = c(1:columns))
    combinations[[1]] <- firstCombi;
    combinations[2:columns] <- lapply(2:columns, function(n) data.frame(t(combn(c(1:columns), n))));  # applies combn for 2:cols
  } else {
    combinations[start:columns] <- lapply(start:columns, function(n) data.frame(t(combn(c(1:columns), n))));
  }

  combinationTable <- do.call(plyr::rbind.fill, combinations);  # combines dataframes by row and fills empty entries with NA
  listOfColumns <- list(scoresTable[,1]);  # add year column
  
  # creates a powerset column for each combination
  combinationNames <- list("year");
  combinations <- list();

  for (i in 1:nrow(combinationTable)) {  # iterate over each combination
    combination <- as.numeric(combinationTable[i,]);
    combination <- combination[!is.na(combination)];  # remove NAs 
    combinations <- rlist::list.append(combinations, combination);
    combinationNames <- rlist::list.append(combinationNames, toString(combination));
    # implements:
    # for (j in 1:nrow(scoresTable)) {  # iterate over each row under this combination
    #   rowSum <- sum(as.numeric(scoresTable[j, combination + 1]));  # +1, since the first column is the year
    #   columnValues <- rlist::list.append(columnValues, rowSum);
    # }
    columnValues <- lapply(1:nrow(scoresTable), function(n) sum(as.numeric(scoresTable[n, combination + 1])));

    listOfColumns <- rlist::list.append(listOfColumns, columnValues);
  }
  
  names(listOfColumns) <- unlist(combinationNames);

  # create data frame out of listOfColumns
  powersetTable <- t(do.call(rbind.data.frame, listOfColumns));
  rownames(powersetTable) <- 1:nrow(powersetTable);
  return(list(scoresTable = as.data.frame(powersetTable), combinations = combinations));
}


#########################################################
#' Returns the best years per column weighted with weighting i.e.
#' weighting often the year is inserted into the returned list,
#' at least the first year - the other can be weighted down.
#'
#' @param bestYearsPerColumn {list} the list of best years to return
#' @param weighting {numerical} tells how to weight the years
#' @param bestYears {vector} the selected years in the column which should
#' @param doubleWeighting {logical} tells if the years with ranks > 1 should be weighted down
#'
#' @return {list} the best years for each column
Analyzer.__getWeightedYearsPerColumn <- function(bestYearsPerColumn, weighting, 
                                                 bestYears, doubleWeighting) {
  if (doubleWeighting) {
    # iterate over all best years
    for (k in 1:length(bestYears)) {
      iterator <- k;  # since k is needed to access the year
      
      if (iterator > weighting) {
        iterator <- weighting;
      }
      
      # iterate over all weightings
      for (j in iterator:weighting) {
        bestYearsPerColumn <- rlist::list.append(bestYearsPerColumn, bestYears[[k]]);
      }
    }
  } else {  # weight every year equally 
    for (j in 1:weighting) {
      bestYearsPerColumn <- rlist::list.append(bestYearsPerColumn, bestYears);
    }
  }
  
  return(bestYearsPerColumn);
}


#########################################################
#' Returns the score for a given rank.
#'
#' @param scoresPerSample {list} the scores per sample
#' @param rank {numeric} the rank for which the score should be extracted
#'
#' @return {list} scores per ranks
#' @export
Analyzer.getRankScore <- function(scoresPerSample, rank) {
  scoresPerRank <- list();
  
  # iterate over each sample
  for (i in 1:length(scoresPerSample)) {
    sampleScores <- scoresPerSample[[i]];
    sortedScores <- sort(unlist(sampleScores));
    scoresPerRank <- rlist::list.append(scoresPerRank, sortedScores[[rank]]);
  }
  
  return(scoresPerRank);
}


#########################################################
#' Compute ranks for a histogram.
#'
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param topYearsPerSample {list} the top years of samples
#'
#' @return {list} ranks per sample
#' @export
Analyzer.getPeakRanks <- function(sampleData, topYearsPerSample) {
  ranks <- list();
  
  # iterate over all samples
  for (i in 1:length(sampleData)) {
    # retrieve
    sample <- sampleData[[i]];
    correctStartYear <- as.numeric(sample$years[[1]]);
    
    topYears <- topYearsPerSample[[i]];
    
    # get histogram for top years
    allYears <- unlist(topYears);
    counts <- table(allYears);

    positionInUniqueYears <- which(allYears == correctStartYear);  # test for existence of correct year in histogram
    
    if (identical(integer(0), positionInUniqueYears)) {  # if not found
      ranks <- rlist::list.append(ranks, Defaults.OUTLIER);
    } else {
      countCorrectYear <- counts[[toString(correctStartYear)]];  # lookup count for correct year in descending order
      sortedCounts <- sort(as.vector(counts), decreasing = TRUE);  # sort the counts in descending order
      ranksWithThatCount <- which(sortedCounts == countCorrectYear);  # lookup positions for correct count in all sorted counts
      peakRank <- ranksWithThatCount[[length(ranksWithThatCount)]];
      ranks <- rlist::list.append(ranks, peakRank);
    }
  }
  
  return(ranks);
}


#########################################################
#' Computes correlation coefficient 
#' at each position of the sample in the chronology.
#'
#' @param consensus {list} consensus (or buckets)
#' @param samples {list} data of multiple curve-files (widths with its years)
#' @param method {string} the correlation method which should be used {"spearman", "pearson", "kendall"}
#' @param minDistanceApproach {logical} tells if the bucket based approach should be applied on the ring-width approach
#' @param inSampleWidths {logical} tells if samples are containers for widths (alternatively the samples are already the widths)
#' 
#' @return {list} the scores per sample
#' @export
Analyzer.computeCoefficientsPerSample <- function(consensus, samples, method, minDistanceApproach, inSampleWidths = FALSE) {
  scoresPerSample <- list();

  if (minDistanceApproach) {
    consensus <- Analyzer.__dataPreparation2(consensus, samples, FALSE, FALSE)$bucketChronology;
  }
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    if (minDistanceApproach && !inSampleWidths) {
      widths <- samples[[i]];
    } else {
      widths <- samples[[i]]$widths;
      
      if(is.null(widths)) {  # extension for maximum densities
        widths <-  samples[[i]]$characteristics;
      }
    }

    scores <- Analyzer.__computeCoefficients(consensus, unlist(widths), method, minDistanceApproach);
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the consensus, which consists out of profiles.
#'
#' @param consensus {list} consensus (widths with its years)
#' @param pattern {vector} the pattern which have to be found in the consensus
#' @param method {string} the correlation method which should be used {"spearman", "pearson", "kendall"}
#' @param minDistanceApproach {logical} tells if the bucket based approach should be applied on the ring-width approach
#' 
#' @return {list} the scores
Analyzer.__computeCoefficients <- function(consensus, pattern, method, minDistanceApproach) {
  scores <- list();
  
  n <- length(consensus);
  m <- length(pattern);
  
  # iterate over each curve position with the subpattern
  for (i in 1:(n-m+1)) {  # "+1", because i is "1" at the beginning and else we wouldn't align against every position 
    # grab a consensus-cut of pattern-size
    cut <- CurvesMiner.getSubset(consensus, i, (i-1) + m);  # "-1", because i is "1" at the beginning
    
    if (minDistanceApproach) {
      # iterate over each position of the pattern and the bucket-cut
      intermediateScores <- list();
      
      # iterate over pattern
      for (k in 1:m) {
        bucket <- cut[[k]];
        widthPattern <- pattern[[k]];
        
        bucketScores <- list();
        
        # iterate over all profiles in the bucket
        for (l in 1:length(bucket)) {
          widthBucket <- bucket[[l]];
          bucketScore <- abs(widthBucket - widthPattern);
          bucketScores <- rlist::list.append(bucketScores, bucketScore);
        }
        
        intermediateScore <- min(unlist(bucketScores));
        intermediateScores <- rlist::list.append(intermediateScores, intermediateScore);
      }
      
      score <- sum(unlist(intermediateScores));
    } else {  # correlation coefficient approach
      
      if (identical(method, Defaults.CORRELATION_GLEICHLAUF)) {
        score <- Math.gleichlauf(unlist(cut), pattern);
      } else if (identical(method, Defaults.CORRELATION_T_VALUE)) {
        score <- Math.twoSampleTtestValues(unlist(cut), pattern);
      } else {
        score <- cor(unlist(cut), pattern, method = method);
      }
    }
    
    scores <- rlist::list.append(scores, score); 
  }
  
  return(scores);
}


#########################################################
#' Returns the scores from the given top years.
#'
#' @param topYears {list|vector} the years from which the scores are from interest
#' @param sampleYears {list} the years from the sample
#' @param sampleScores {list} the scores from the sample
#'
#' @return {list} the score from the top years
Analyzer.getScoresFromYears <- function(topYears, sampleYears, sampleScores) {
  scores <- list();
  
  # iterate over all years and look up scores
  for (i in 1:length(topYears)) {
    year <- topYears[[i]];
    position <- which(sampleYears == year);  
    correctScore <- sampleScores[[position]];
    
    scores <- rlist::list.append(scores, correctScore);
  }
  
  return(scores);
}


#########################################################
#' Computes the top rated samples in the ring-width based approach.
#'
#' @param number {numerical} the number of samples on which you want to look on
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} contains different scoreTypes
#'
#' @return {list(fileNames, ranks)} the filenames of top rated samples together with their ranks
#' @export
Analyzer.getTopSamples <- function(number, passPathGenericName, 
                                   scoresFolderGenericName, 
                                   passes, scoreType) {
  # compute ranks out of scores
  ranksPerPass <- list();
  sampleNamesPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    ranksPerScoreType <- list();
    
    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # read in score-type scores per sample
    scoresData <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                            scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR, TRUE);
    
    # get scores for correct position of the samples
    scoresPerSample <- scoresData$scoresData;
    scoresForCorrectPositions <- Analyzer.getScoresForCorrectSamplePositions(scoresPerSample, sampleData,
                                                                             masterChronology$years);
    
    # for each sample: compute the rank of the correct-position score within all scores of the sample
    ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, TRUE);
    names <- scoresData$fileNames;
    
    # store the ranks per pass
    ranksPerPass <- rlist::list.append(ranksPerPass, ranks);
    sampleNamesPerPass <- rlist::list.append(sampleNamesPerPass, names);
  }
  
  # unlist both get order of ranks
  ranks <- unlist(ranksPerPass);
  names <- unlist(sampleNamesPerPass);
  orderRanks <- order(ranks);
  
  sortedRanks <- ranks[orderRanks];
  names <- names[orderRanks];
  
  return(list(fileNames = names[1:number], ranks = sortedRanks[1:number]));
}


#########################################################
#' Given filepaths, the bucket based approach is applied 
#' and the ranks are looked up and returned.
#'
#' @param filepaths {list} the paths to the files
#'
#' @return {list} the ranks from the files that have been loaded
#' @export
Analyzer.getFilesRanks <- function(filepaths) {
  ranks <- list();
  
  # iterate over all paths
  for (i in 1:length(filepaths)) {
    # load data
    filepath <- filepaths[[i]];
    scoreData <- read.csv(filepath, header = TRUE, sep = Symbols.REAL_DATA_SEPARATOR);
    correctYears <- stringr::str_extract(filepath, Expression.YEARS_IN_BETWEEN);
    firstYear <- as.numeric(stringr::str_split(correctYears, Symbols.HYPHEN)[[1]][[1]]);
    
    # lookup correct score within scores
    scores <- scoreData$score;
    years <- scoreData$year;
    
    correctYearIndex <- which(years == firstYear);
    correctScore <- scores[[correctYearIndex]];
    
    # sort the scores and compute the rank of the correct score
    sortedScores <- sort(scores);
    ranksWithThatScore <- which(sortedScores == correctScore);
    rankOfCorrectScore <- ranksWithThatScore[[length(ranksWithThatScore)]];  # last one
    
    # store the rank
    ranks <- rlist::list.append(ranks, rankOfCorrectScore);
  }
  
  return(ranks);
}


#########################################################
#' Computes the buckets-based corresponding ranks for ring-width-approach.
#'
#' @param topYears {numerical} the number of top scores that have to be retrieved 
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the score-type for the ring-width approach
#' @param histogramScoreType {list} the scoreType of the buckets-based method
#'
#' @return {list(outcasts = outcasts, ranks = ranks)} outcasts and ranks
#' @export
Analyzer.getCorrespondingRanksAndOutcasts <- function(topYears, passPathGenericName, 
                                                      scoresFolderGenericName, passes, 
                                                      scoreType, histogramScoreType) {
  outcasts <- 0;
  ranks <- list();
  
  # the top years per pass and file 
  dataPerPass <- Analyzer.getTopYearsPerPass(topYears, passPathGenericName, 
                                             scoresFolderGenericName, passes, 
                                             scoreType);
  
  topYearsPerPass <- dataPerPass$topYearsPerPass;
  filePathsPerPass <- dataPerPass$filePathsPerPass;
  
  # iterate over each pass
  for (i in 1:length(topYearsPerPass)) {
    passTopYears <- topYearsPerPass[[i]];
    passFilePaths <- filePathsPerPass[[i]];
    
    # convert names to load corresponding scores in buckets-based approach
    finalReplacement4 <- paste(Strings.SCORE, Symbols.UNDERSCORE, histogramScoreType, sep = Symbols.EMPTY);
    finalReplacement5 <- paste(Strings.SCORES, Symbols.UNDERSCORE, histogramScoreType, sep = Symbols.EMPTY);
    
    passFilePaths <- gsub(Strings.TO_REPLACE_3, Strings.REPLACEMENT_3, passFilePaths);
    passFilePaths <- gsub(Strings.TO_REPLACE_6, Strings.REPLACEMENT_3, passFilePaths);
    passFilePaths <- gsub(Strings.TO_REPLACE_4, finalReplacement5, passFilePaths);
    passFilePaths <- gsub(Strings.TO_REPLACE_5, finalReplacement4, passFilePaths);
    passFilePaths <- gsub(Strings.TO_REPLACE_7, finalReplacement5, passFilePaths);
    passFilePaths <- gsub(Strings.TO_REPLACE_8, finalReplacement4, passFilePaths);
    
    # iterate over each sample
    for (j in 1:length(passFilePaths)) {
      filePath <- passFilePaths[[j]];
      
      # load file
      table <- read.csv(filePath, header = TRUE, sep = Symbols.REAL_DATA_SEPARATOR);
      years <- table$year;
      scores <- table$score;
      
      # look up correct score (from the given startyear written on file)
      correctYears <- stringr::str_extract(filePath, Expression.YEARS_IN_BETWEEN);
      correctStartYear <- as.numeric(stringr::str_split(correctYears, Symbols.HYPHEN)[[1]][[1]]);
      position <- which(years == correctStartYear);  
      correctScore <- scores[[position]];
      
      # look up top years scores in buckets-approach
      topYearsScores <- Analyzer.getScoresFromYears(unlist(passTopYears[[j]]), years, scores);
      sortedTopYearScores <- sort(unlist(topYearsScores));
      
      # look which rank that correct score has in the upper sorted scores
      # and if it is included or an outcast
      ranksWithThatScore <- which(sortedTopYearScores == correctScore);  # position of correct score
      
      if (identical(integer(0), ranksWithThatScore)) {
        outcasts <- outcasts + 1;
      } else {
        rankOfCorrectScore <- ranksWithThatScore[[length(ranksWithThatScore)]];  # last one
        ranks <- rlist::list.append(ranks, rankOfCorrectScore);
      }
    }
  }
  
  return(list(outcasts = outcasts, ranks = ranks));
}


#########################################################
#' Return the top years per sample from several passes.
#'
#' @param topYears {numerical} the number of top scores that have to be retrieved 
#' for each sample with the ring-width based approach
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {list} contains different scoreTypes
#'
#' @return {list(topYearsPerPass, filePathsPerPass)} the top years per sample together with
#' @export
Analyzer.getTopYearsPerPass <- function(topYears, passPathGenericName, 
                                        scoresFolderGenericName, passes, 
                                        scoreType) {
  # get top years e.g. top 20 from ringwidth-based approach for each sample
  # to limit the number of years that have to be checked with the buckets-based approach
  topYearsPerPass <- list();
  filePathsPerPass <- list();
  
  # iterate over each pass to find the top years for each sample
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);

    sampleData <- Loader.readInCurves(paste(passesPath, Paths.OSTALB_TESTSET_FOLDER, sep = Paths.PATH), 
                                      Symbols.REAL_DATA_SEPARATOR);
    
    masterChronology <- Loader.readInProfiles(paste(passesPath, Paths.PATH ,sep = Symbols.EMPTY), 
                                              Files.OSTALB_CHRONOLOGY, Symbols.REAL_DATA_SEPARATOR);
    
    # read in score-type scores per sample
    scoresData <- Loader.readInScores(paste(passesPath, Paths.PATH, scoresFolderGenericName, 
                                            scoreType, sep = Symbols.EMPTY), Symbols.REAL_DATA_SEPARATOR, TRUE);
    
    # store top rated positions or years with the highest Pearson Correlation Coefficients
    topYearsSample <- Analyzer.getTopYearsPerSample(scoresData$scoresTables, topYears, FALSE, 
                                                    FALSE, FALSE, TRUE);
    
    # store the ranks per pass
    topYearsPerPass <- rlist::list.append(topYearsPerPass, topYearsSample);
    filePathsPerPass <- rlist::list.append(filePathsPerPass, scoresData$fileNames);
  }
  
  return(list(topYearsPerPass = topYearsPerPass, filePathsPerPass = filePathsPerPass));
}


#########################################################
#' Returns the scores per sample from several passes.
#' Hint: This function assumes that all scores are unique.
#'
#' @param scoresPerSample {list} per sample scores
#' @param scoresToChange {list} the top score values
#' @param weights {vector} the delta scores
#'
#' @return {list} scoresPerSample with changed years
Analyzer.changeSampleScore <- function(scoresPerSample, scoresToChange, weights) {
  changedScoresPerSample <- list();
  
  # iterate over each list of scores
  for (i in 1:length(scoresPerSample)) {
    print(i);
    sampleScores <- scoresPerSample[[i]];
    scoreToChange <- scoresToChange[[i]];
    weight <- weights[[i]];
    
    if (weight < 1) {
      # look up where the score is that should be changed
      position <- which(sampleScores == scoreToChange);
      sampleScores[[position]] <- sampleScores[[position]] * 1/weight;
    }
    
    changedScoresPerSample <- rlist::list.append(changedScoresPerSample, sampleScores);
  }
  
  return(changedScoresPerSample);
}


#########################################################
#' Returns the top-years per sample.
#'
#' @param scoresPerSample {list} the scores for each sample
#' @param sampleData {list} data of multiple curve-files (profiles and widths with its years)
#' @param years {vector} the years of the scores (each score has a corresponding year)
#' @param bestYearsMax {numerical} tells how many top years per sample should be returned
#' @param qualities {list} tells which quality measures should be shown in the encoded matrix
#'
#' @return {matrix} the years per sample
#' @example 
#' #matrix
#' sample,  rank[1],  rank[2], ..., rank[n]
#' 401,     2015,     2017,     .., 1922
#' 402,     2013,     2030,     .., 1629
#' ...
#' 
#' @export
Analyzer.getYearsPerSample <- function(scoresPerSample, sampleData, years, 
                                       bestYearsMax, qualities) {
  yearsPerSample <- list();
  pValuesPerSample <- list();
  distancesPerSample <- list();
  
  names <- list();
  
  # iterate over each sample
  for (i in 1:length(scoresPerSample)) {
    scoresOfSample <- unlist(scoresPerSample[[i]]);
    name <- sampleData[[i]]$name;
      
    orderOfScores <- order(scoresOfSample);
    
    if (length(qualities) > 0) {
      if (qualities$pValues) {  # p-values per sample
        pValue <- Analyzer.getPValue(scoresOfSample, 1);
        pValuesPerSample <- rlist::list.append(pValuesPerSample, pValue);
      } 
      
      if (qualities$scores) {  # scores per sample
        sortedScores <- scoresOfSample[orderOfScores];
        bestScores <- round(sortedScores[1:bestYearsMax], digits = Defaults.ROUNDING_DIGITS);
        distancesPerSample <- rlist::list.append(distancesPerSample, bestScores);
      }
    }

    bestYears <- Analyzer.__getYearsOfSample(orderOfScores, years, bestYearsMax);
    
    yearsPerSample <- rlist::list.append(yearsPerSample, bestYears);
    names <- rlist::list.append(names, name);
  }
  
  return(Analyzer.__getMatrix(yearsPerSample, names, pValuesPerSample, distancesPerSample));
}


#########################################################
#' Returns the top-years of a sample with corresponding quality measures.
#'
#' @param orderOfScores {vector} the sorting for the scores
#' @param years {vector} the years of the scores (each score has a corresponding year)
#' @param bestYearsMax {numerical} tells how many top years per sample should be returned
#'
#' @return {list} the top years of a sample
#' @export
Analyzer.__getYearsOfSample <- function(orderOfScores, years, bestYearsMax) {
  sortedYears <- years[orderOfScores];
  
  numYears <- length(sortedYears);
  
  if (bestYearsMax <= numYears && bestYearsMax > 0) {
    return(sortedYears[1:bestYearsMax]);
  }
  
  return(sortedYears[1:numYears]);
}


#########################################################
#' Returns a matrix which shows per sample the most promising years.
#'
#' @param yearsPerSample {list} the sorted years per sample
#' @param names {list} the names of the sample
#' @param pValuesPerSample {list} the p-values per sample 
#' @param scoresPerSample {list} the best scores per sample
#'
#' @return {matrix} the name matrix
#' @export
Analyzer.__getMatrix <- function(yearsPerSample, names, 
                                 pValuesPerSample, scoresPerSample) {
  
  addPvalues <- length(pValuesPerSample) != 0;
  addScores <- length(scoresPerSample) != 0;

  # combine yearsPerSample per row
  if (addScores) {
    perSampleYears <- do.call(rbind, yearsPerSample);
    perSampleScores <- do.call(rbind, scoresPerSample);
    
    columns <- list();
    
    for (i in 1:ncol(perSampleYears)) {  # place score right of rank
      years <- perSampleYears[,i];
      scores <- perSampleScores[,i];
      
      # append
      columns <- rlist::list.append(columns, years);
      columns <- rlist::list.append(columns, scores);
    }
    
    perSampleYears <- do.call(cbind, columns);
  } else {
    perSampleYears <- do.call(rbind, yearsPerSample);
  }

  # add p-value column if necessary
  if (addPvalues) {
    values <- round(unlist(pValuesPerSample), digits = Defaults.ROUNDING_DIGITS_2);
    perSampleYears <- cbind(values, perSampleYears);
  }
  
  # add column with names
  perSampleYears <- cbind(unlist(names), perSampleYears);
  
  colnames(perSampleYears) <- Analyzer.__getColNames(length(yearsPerSample[[1]]), addPvalues, addScores);
  
  return(perSampleYears);
}


#########################################################
#' Returns the top-years of a sample.
#'
#' @param ranks {numerical} the number of ranks to display
#' @param addPvalues {logical} tells if p-value columns should be added or not
#' @param addScores {logical} tells if scores columns should be added or not
#'
#' @return {vector} the colnames
#' @export
Analyzer.__getColNames <- function(ranks, addPvalues, addScores) {
  colnames <- list(Defaults.SAMPLE_NAME);
  
  if (addPvalues) {
    colnames <- rlist::list.append(colnames, Defaults.P_VALUE_NAME);
  }
  
  for (i in 1:ranks) {
    name <- paste(Defaults.RANK_NAME, i, sep = Symbols.EMPTY);
    colnames <- rlist::list.append(colnames, name);
    
    if (addScores) {
      name <- paste(Defaults.SCORE_NAME, i, sep = Symbols.EMPTY);
      colnames <- rlist::list.append(colnames, name);
    }
  }
  
  return(unlist(colnames));
}


#########################################################
#' Creates scores tables out of the per sample scores.
#'
#' @param scoresPerSample {list} the scores per sample
#' @param years {list} the years from which that scores
#'
#' @return {list} the scores tables
#' @export
Analyzer.getScoresTables <- function(scoresPerSample, years) {
  tables <- list();
  
  for (i in 1:length(scoresPerSample)) {
    # retrieve scores and sample-name
    scores <- scoresPerSample[[i]];
    
    # create file with that sample-name
    data <- do.call(rbind, scores);
    data <- cbind(years[1:nrow(data)], data);
    
    names <- list();
    
    # if one score, then write just "score", else "score[i]" i.e. "score1" and so on
    if (ncol(data) != 2) {
      for (j in 1:(ncol(data) - 1)) {
        name <- paste(Strings.SCORE_LOWER_CASE, j, sep = Symbols.EMPTY);
        names <- rlist::list.append(names, name);
      }
    } else {
      name <- Strings.SCORE_LOWER_CASE;
      names <- rlist::list.append(names, name);
    }
    
    colnames(data) <- c(Strings.YEAR_LOWER_CASE, unlist(names));
    
    tables <- rlist::list.append(tables, as.data.frame(data));
  }
  
  return(tables);
}


#########################################################
#' Returns the top-years per sample (for voting approach).
#'
#' @param topYearsPerSample {list} the best years per sample
#' @param sampleData {list} data of multiple curve-files (profiles and widths with its years)
#' @param bestYearsMax {numerical} tells how many years should be stored
#'
#' @return {matrix} the years per sample
#' @example 
#' #matrix
#' sample,  rank[1],  rank[2], ..., rank[n]
#' 401,     2015,     2017,     .., 1922
#' 402,     2013,     2030,     .., 1629
#' ...
#' 
#' @export
Analyzer.getYearsPerSample2 <- function(topYearsPerSample, sampleData, bestYearsMax) {
  yearsPerSample <- list();
  names <- list();
  
  # iterate over each sample
  for (i in 1:length(topYearsPerSample)) {
    topYears <- topYearsPerSample[[i]];
    name <- sampleData[[i]]$name;
    
    yearsOfSample <- Analyzer.__getYearsOfSample2(unlist(topYears), bestYearsMax);
    
    # fill with NAs if not enough best years available
    numOfNA <- bestYearsMax - length(yearsOfSample);
    yearsOfSample <- c(yearsOfSample, rep(NA, numOfNA));
    
    yearsPerSample <- rlist::list.append(yearsPerSample, yearsOfSample);
    names <- rlist::list.append(names, name);
  }
  
  # create matrix
  # combine yearsPerSample per row
  perSampleYears <- do.call(rbind, yearsPerSample);
  
  # add column with names
  perSampleYears <- cbind(unlist(names), perSampleYears);
  
  colnames(perSampleYears) <- Analyzer.__getColNames(bestYearsMax, FALSE, FALSE);
  
  return(perSampleYears);
}


#########################################################
#' Returns the top-years of a sample (for voting approach).
#'
#' @param topYearsOfSample {vector} the best years per sample
#' @param bestYearsMax {numerical} tells how many years should be stored
#'
#' @return {list} the top years of a sample
#' @export
Analyzer.__getYearsOfSample2 <- function(topYearsOfSample, bestYearsMax) {
  # create histogram
  years <- as.numeric(levels(factor(topYearsOfSample)));  # years
  counts <- as.vector(table(topYearsOfSample));  # corresponding count

  orderOfCounts <- order(counts, decreasing = TRUE);  # order of counts -> sorting histogram
  sortedYears <- years[orderOfCounts];
  
  # return data
  numYears <- length(sortedYears);
  
  if (bestYearsMax <= numYears && bestYearsMax > 0) {
    return(sortedYears[1:bestYearsMax]);
  }
  
  return(sortedYears[1:numYears]);
}


#########################################################
#' Returns a table which returns for a year the index within all years.
#'
#' @param years {vector} the years for which a table has to be created
#'
#' @return {environment} indices per year structure
#' @export
Analyzer.getIndicesPerYear <- function(years) {
  indicesPerYear <- new.env();
  
  i <- 1;
  
  for (y in years) {
    indicesPerYear[[toString(y)]] <- i;
    i <- i + 1;
  }
  
  return(indicesPerYear);
}


#########################################################
#' Returns the top-years per sample.
#'
#' @param scoresPerSample {list} the scores for each sample
#' @param sampleData {list} data of multiple curve-files (profiles and widths with its years)
#' @param yearsPerSample {vector} the years of the scores (each score has a corresponding year)
#' @param bestYearsMax {numerical} tells how many top years per sample should be returned
#' @param qualities {list} tells which quality measures should be shown in the encoded matrix
#'
#' @return {matrix} the years per sample
#' @example 
#' #matrix
#' sample,  rank[1],  rank[2], ..., rank[n]
#' 401,     2015,     2017,     .., 1922
#' 402,     2013,     2030,     .., 1629
#' ...
#' 
#' @return {list} the years per sample
#' @export
Analyzer.getYearsPerSample3 <- function(scoresPerSample, sampleData, 
                                        yearsPerSample, bestYearsMax, qualities) {
  sortedYearsPerSample <- list();
  names <- list();
  pValuesPerSample <- list();
  distancesPerSample <- list();
  
  # iterate over each sample
  for (i in 1:length(scoresPerSample)) {
    scoresOfSample <- unlist(scoresPerSample[[i]]);
    orderOfScores <- order(scoresOfSample);
    yearsOfSample <- yearsPerSample[[i]];
      
    name <- sampleData[[i]]$name;
    
    if (qualities$pValues) {  # p-values per sample
      pValue <- Analyzer.getPValue(scoresOfSample, 1);
      pValuesPerSample <- rlist::list.append(pValuesPerSample, pValue);
    } 
    
    if (qualities$scores) {  # scores per sample
      sortedScores <- scoresOfSample[orderOfScores];
      bestScores <- round(sortedScores[1:bestYearsMax], digits = Defaults.ROUNDING_DIGITS);
      distancesPerSample <- rlist::list.append(distancesPerSample, bestScores);
    }
    
    sortedYears <- Analyzer.__getYearsOfSample(orderOfScores, unlist(yearsOfSample), bestYearsMax);
    
    sortedYearsPerSample <- rlist::list.append(sortedYearsPerSample, sortedYears);
    names <- rlist::list.append(names, name);
  }
  
  return(Analyzer.__getMatrix(sortedYearsPerSample, names, 
                              pValuesPerSample, distancesPerSample));
}


#########################################################
#' Moves through the chronology and computes the correlation coefficient
#' between neighboured consensus-profiles. All such computed scores
#' are stored and finally returned
#'
#' @param masterChronology {list} the list of profiles (can have different amount of points)
#'
#' @return {list} scores between neighboured profiles which were interpolated to same length
#' @export
Analyzer.compareNeighbouredProfiles <- function(masterChronology) {
  n <- length(masterChronology);
  
  neighbourScores <- list();
  
  # iterate over all profiles
  for (i in 1:(n-1)) {
    profileLeft <- masterChronology[[i]];
    profileRight <- masterChronology[[i+1]];
    
    # interpolation of x-values is necessary i.e. both curves should have the same number
    # Hint: does exactly the same like MICA interpolation function, but significantly faster
    profileLeft <- unlist(Math.interpolatePoints(profileLeft, Parameters.RANKING_DIST_SAMPLE));
    profileRight <- unlist(Math.interpolatePoints(profileRight, Parameters.RANKING_DIST_SAMPLE));
    
    score <- cor(profileLeft, profileRight);
    neighbourScores <- rlist::list.append(neighbourScores, score);
  }
  
  return(neighbourScores);
}


#########################################################
#' Return the potential pointer-years or event-years by retrieving all
#' years with scores lower the limit.
#'
#' @param scores {list} the correlation coefficients
#' @param years {vector} the corresponding years for the scores
#' @param limit {numerical} the limit at which all scores below identified as characteristic
#' @param indices {logical} tells if instead of years the corresponding indices within the sample
#' had to be returned
#'
#' @return {list} the characteristic years
Analyzer.getCharacteristicYears <- function(scores, years, limit, indices = FALSE) {
  characteristicYears <- list();
  
  # iterate over all scores
  for (i in 1:length(scores)) {
    score <- scores[[i]];
    
    if (score < limit) {
      if (indices) {
        characteristicYears <- rlist::list.append(characteristicYears, i + 1);  # +1, since the first year cannot be detected as event-year 
      } else {
        characteristicYears <- rlist::list.append(characteristicYears, years[[i]]);
      }
    }
  }
  
  return(characteristicYears);
} 


#########################################################
#' Divides the sample data into characteristic and non-characteristic years.
#'
#' @param sampleData {list} data of multiple curve-files (profiles and widths with its years)
#' @param correlationCoefficientLimit {numerical} the limit below which all years defined as pointer-years
#'
#' @return {list(samplesWithPointerYears, characteristicYearsPerSample, samplesWithoutPointerYears)} the samples split
#' @export
Analyzer.splitIntoCharacteristicYears <- function(sampleData, correlationCoefficientLimit) {
  samplesWithCharacteristicYears <- list();
  samplesWithoutCharacteristicYears <- list();
  
  characteristicYearsIndicesPerSample <- list();
  
  # iterate over each sample
  for (i in 1:length(sampleData)) {
    sample <- sampleData[[i]];
    
    neighbourScores <- Analyzer.compareNeighbouredProfiles(sample$profiles);
    
    # first year if characteristic, not detectable
    characteristicStartYears <- Analyzer.getCharacteristicYears(neighbourScores, 
                                                           unlist(sample$years[2:length(sample$years)]),  
                                                           correlationCoefficientLimit, TRUE);
    # if characteristic years contained, store
    if (length(characteristicStartYears) > 0) {
      samplesWithCharacteristicYears <- rlist::list.append(samplesWithCharacteristicYears, sample);
      characteristicYearsIndicesPerSample <- rlist::list.append(characteristicYearsIndicesPerSample, characteristicStartYears);
    } else {
      samplesWithoutCharacteristicYears <- rlist::list.append(samplesWithoutCharacteristicYears, sample);
    }
  }
  
  return(list(samplesWithCharacteristicYears = samplesWithCharacteristicYears,
              characteristicYearsIndicesPerSample = characteristicYearsIndicesPerSample,
              samplesWithoutCharacteristicYears = samplesWithoutCharacteristicYears));
}


#########################################################
#' Computes the ranks for samples without event years.
#'
#' @param allYears {vector} the vector of all years
#' @param scoresPath {string} the path to the samples scores
#' @param samplesWithoutEventYears {list} samples without event years
#' @param scoreType {list} the scoreType of the ring-width approach
#'
#' @return {list} ranks of the given samples
#' @export
Analyzer.getRanksForUsualSamples <- function(allYears, scoresPath, 
                                             samplesWithoutEventYears, 
                                             scoreType) {
  allRanks <- list();
  
  # iterate over each sample without an event-year
  for (j in 1:length(samplesWithoutEventYears)) { 
    # load scores given the score type and sample name from above
    sample <- samplesWithoutEventYears[[j]];
    name <- sample$name;
    scoresToLoad <- paste(scoresPath, name, Strings.SCORE, 
                          Symbols.UNDERSCORE, scoreType, Extensions.CSV, sep = Symbols.EMPTY);
    
    # read in score-type scores per sample
    entries <- read.csv(scoresToLoad, header = TRUE, sep = Symbols.REAL_DATA_SEPARATOR);
    scores <- entries$score;
    
    # look up correct score for the sample
    correctStartYear <- as.numeric(sample$years[[1]]);
    index <- which(allYears == correctStartYear);
    correctScore <- scores[[index]];
    
    # look which rank that correct score has
    sortedScores <- sort(scores);
    ranksWithThatScore <- which(sortedScores == correctScore);  
    finalRank <- ranksWithThatScore[[length(ranksWithThatScore)]];  # not necessary usually, since scores should be unqiue
    allRanks <- rlist::list.append(allRanks, finalRank);
  }
  
  return(allRanks)
}


#########################################################
#' Computes the correct years ranks for samples with event years.
#'
#' @param allYears {vector} the vector of all years
#' @param scoresPath {string} the path to the samples scores
#' @param samplesWithEventYears {list} samples with event years
#' @param eventYearIndicesPerSample {list} the indices which tell 
#' where in the samples is the event-year
#' @param scoreType {list} the scoreType of the ring-width approach
#' @param pointerYears {list} the pointer-years from an chronology
#' 
#' @return {list(predictedRanks, outcasts)} ranks of the given samples correct year
#' @export
Analyzer.getRanksForCharacteristicSamples <- function(allYears, scoresPath, samplesWithEventYears, eventYearIndicesPerSample,
                                                      scoreType, pointerYears) {
  
  predictedScoresPerSample <- Analyzer.__getPredictedScores(scoresPath, samplesWithEventYears, eventYearIndicesPerSample,
                                                            scoreType, pointerYears);
  predictedRanks <- list();
  outcasts <- 0;
  
  # iterate over all samples
  for (i in 1:length(predictedScoresPerSample)) {
    scoresOfSample <- unlist(predictedScoresPerSample[[i]]);
    sortedScoresOfSample <- sort(scoresOfSample);

    # load scores given the score type and sample name from above
    sample <- samplesWithEventYears[[i]];
    name <- sample$name;
    scoresToLoad <- paste(scoresPath, name, Strings.SCORE, 
                          Symbols.UNDERSCORE, scoreType, Extensions.CSV, sep = Symbols.EMPTY);
    
    # read in score-type scores per sample
    entries <- read.csv(scoresToLoad, header = TRUE, sep = Symbols.REAL_DATA_SEPARATOR);
    scores <- entries$score;
    
    # look up rank for correct score
    correctStartYear <- as.numeric(sample$years[[1]]);
    index <- which(allYears == correctStartYear);
    correctScore <- scores[[index]];
    
    # look which rank that correct score has or if it is an outlier
    ranksWithThatScore <- which(sortedScoresOfSample == correctScore); 
    
    if (identical(ranksWithThatScore, integer(0))) {  # if correctScore not contained
      outcasts <- outcasts + 1;
    } else {
      finalRank <- ranksWithThatScore[[length(ranksWithThatScore)]];
      predictedRanks <- rlist::list.append(predictedRanks, finalRank);
    }
  }
  
  return(list(predictedRanks = predictedRanks, 
              outcasts = outcasts));
}


#########################################################
#' Computes the predicted scores for samples with event years.
#'
#' @param scoresPath {string} the path to the samples scores
#' @param samplesWithEventYears {list} samples with event years
#' @param eventYearIndicesPerSample {list} the indices which tell 
#' where in the samples is the event-year
#' @param scoreType {list} the scoreType of the ring-width approach
#' @param pointerYears {list} the pointer-years from an chronology
#'
#' @return {list} the predicted scores
#' @export
Analyzer.__getPredictedScores <- function(scoresPath, samplesWithEventYears, eventYearIndicesPerSample,
                                          scoreType, pointerYears) {
  predictedScoresPerSample <- list();
  
  # iterate over each sample with an event-year
  for (j in 1:length(samplesWithEventYears)) {
    # load scores given the score type and sample name from above
    sample <- samplesWithEventYears[[j]];
    
    sampleName <- sample$name;
    scoresToLoad <- paste(scoresPath, sampleName, Strings.SCORE, 
                          Symbols.UNDERSCORE, scoreType, Extensions.CSV, sep = Symbols.EMPTY);
    
    entries <- read.csv(scoresToLoad, header = TRUE, sep = Symbols.REAL_DATA_SEPARATOR);
    scores <- entries$score;
    years <- entries$year;
    
    eventYearIndicesSample <- eventYearIndicesPerSample[[j]];
    startYearsForSample <- Analyzer.__findStartYearForSample(pointerYears, eventYearIndicesSample);
    
    predictedScoresOfSample <- list();
    
    # lookup the scores for start years of the sample
    # iterate over all years
    for (k in 1:length(startYearsForSample)) {
      year <- startYearsForSample[[k]];
      index <- which(years == year);
      if (!identical(index, integer(0))) {  # it can happen that start-year is outside the range of possible years (due to length of sample)
        score <- scores[[index]];
        predictedScoresOfSample <- rlist::list.append(predictedScoresOfSample, score);
      }
    }
    
    predictedScoresPerSample <- rlist::list.append(predictedScoresPerSample, predictedScoresOfSample);
  }
  
  return(predictedScoresPerSample);
}


#########################################################
#' Finds predicted start years for a sample.
#'
#' @param pointerYears {list} the pointer-years from an chronology
#' @param eventYearIndicesSample {list} the indices which tell 
#' where in the sample is the event-year
#'
#' @return {list} the predicted start years (unique)
#' @export
Analyzer.__findStartYearForSample <- function(pointerYears, eventYearIndicesSample) {
  startYearsForSample <- list();
  
  # iterate over all found pointer-years
  for (k in 1:length(pointerYears)) {
    pointerYear <- as.numeric(pointerYears[[k]]);
    
    # iterate over all indices of the sample
    for (l in 1:length(eventYearIndicesSample)) {
      index <- eventYearIndicesSample[[l]];
      startYear <- pointerYear - index + 1;  # +1, i.e. 1992 is pointer-year and position 2 in sample is event-year, then 1991 is startyear
      startYearsForSample <- rlist::list.append(startYearsForSample, startYear);
    }
  }
  
  return(unique(unlist(startYearsForSample)));
}


#########################################################
#' Computes the weighted difference in the votes of the best 
#' and second best year (in the histogram of years).
#' Weighted means, that the number of second best is counted and
#' the peak is divivded by this number.
#'
#' @param topYearsPerSample {list} the top years of samples
#'
#' @return {list} the delta-peaks of first and second peak per sample
#' @example 
#' years:   1945  1992  1993  1996  2012
#' counts:  15    13    13    13    11
#' 
#' then: 
#' diffPeak = 2
#' weight = 3
#' 
#' final diff: 2/3
#' 
#' @return {list} weights of peaks
#' @export
Analyzer.getWeightedDeltaPeaks <- function(topYearsPerSample) {
  deltaPeaks <- list();
  
  # iterate over all samples
  for (i in 1:length(topYearsPerSample)) {
    topYears <- topYearsPerSample[[i]];
    
    # get histogram for top years
    allYears <- unlist(topYears);
    counts <- as.vector(table(allYears));  # histogram
    
    sortedCounts <- sort(counts, decreasing = TRUE);
    
    if (length(sortedCounts) > 1) {
      weight <- which(sortedCounts == sortedCounts[[2]]);
      deltaPeak <- sortedCounts[[1]] - sortedCounts[[2]];
      deltaPeaks <- rlist::list.append(deltaPeaks, deltaPeak/weight);  # if there are multiple with same weight
    } else {
      deltaPeaks <- rlist::list.append(deltaPeaks, Defaults.ALL_VOTED_FOR_THE_SAME);  # -1, all voted for the same
    }
  }
  
  return(deltaPeaks);
}


#########################################################
#' Computes the p-values per sample.
#'
#' @param scoresPerSample {list} all computed scores of a sample
#' @param threshold {numeric} the number of lowest scores to look at
#' @param deltaScores {logical} tells if the right sided delta-score p-value should be computed 
#' 
#' @return {list} the p-values per sample
#' @export
Analyzer.getPValues <- function(scoresPerSample, threshold, computeDeltaScores = FALSE) {
  pValuesPerSample <- list();
  
  # iterate over all samples
  for (i in 1:length(scoresPerSample)) {
    scoresOfSample <- unlist(scoresPerSample[[i]]);
    
    pValue <- -1;
    
    if (computeDeltaScores) {
      pValue <- Analyzer.getRightTailedPValue(scoresOfSample);
    } else{
      pValue <- Analyzer.getPValue(scoresOfSample, threshold);
    }

    pValuesPerSample <- rlist::list.append(pValuesPerSample, pValue);
  }
  
  return(pValuesPerSample);
}


#########################################################
#' Computes the right p-value from a set of scores and a sample.
#'
#' @param scores {vector} the scores from which a distribution is build 
#' @param threshold {numeric} the number of lowest scores to look at
#'
#' @return {numeric} the p-value
#' @export
Analyzer.getRightTailedPValue <- function(scores) {
  # hint: theoretically it can be used the minimum-function instead of sorting-function, but depending on distance function,
  # it could happen that we have two times same score, so this is the safe variant
  sortedScores <- sort(scores);
  bestScore <- sortedScores[[1]];
  remainingScores <- sortedScores[-1];  # delete first one, since it is looked on differences to the best value
  deltaScores <- remainingScores - bestScore;

  fit <- fitdistrplus::fitdist(deltaScores, "gamma");

  # ########################################################################################################
  # # check which is best distribution for data
  # fitW <- fitdistrplus::fitdist(deltaScores, "weibull");
  # fitg <- fitdistrplus::fitdist(deltaScores, "gamma");
  # fitln <- fitdistrplus::fitdist(deltaScores, "lnorm");
  # fitdistrplus::denscomp(list(fitW, fitg, fitln), legendtext = c("Weibull", "gamma", "lognormal"));
  # 
  # # nicer ggplots
  # # weibullPlot <- stat_function(fun = dweibull, args = list(shape = fitW$estimate[1], scale = fitW$estimate[2]),
  # #                              color = Colors.RED);
  # #
  # # gammaPlot <- stat_function(fun = dgamma, args = list(shape = fitg$estimate[1], rate = fitg$estimate[2]),
  # #                            color = Colors.VIOLET);
  # #
  # # logNormalPlotAndHist <- Plotter.getLogNormalDistributionPlot(scores, fitln, Colors.BLUE, print = FALSE);
  # #
  # # plot <- logNormalPlotAndHist + weibullPlot + gammaPlot;
  # #
  # # print(plot);
  # 
  # data <- fitdistrplus::gofstat(list(fitW, fitg, fitln), fitnames = c("weibull", "gamma", "lognormal"))
  # print(as.character(data$kstest["weibull"]) == "not rejected");
  #########################################################################################################

  pValue <- pgamma(deltaScores[[length(deltaScores)]], fit$estimate[1], fit$estimate[2], lower.tail = FALSE);  # P[X > threshold]
  return(pValue);
}


#########################################################
#' Computes the p-value from a set of scores and a sample.
#'
#' @param scores {vector} the scores from which a distribution is build 
#' @param threshold {numeric} the number of lowest scores to look at
#'
#' @return {numeric} the p-value
#' @export
Analyzer.getPValue <- function(scores, threshold) {
  # # check which is best distribution for data
  # fitW <- fitdistrplus::fitdist(scores, "weibull");
  # fitg <- fitdistrplus::fitdist(scores, "gamma");
  # fitln <- fitdistrplus::fitdist(scores, "lnorm");
  # fitdistrplus::denscomp(list(fitW, fitg, fitln), legendtext = c("Weibull", "gamma", "lognormal"));
  # 
  # # nicer ggplots
  # # weibullPlot <- stat_function(fun = dweibull, args = list(shape = fitW$estimate[1], scale = fitW$estimate[2]),
  # #                              color = Colors.RED);
  # #
  # # gammaPlot <- stat_function(fun = dgamma, args = list(shape = fitg$estimate[1], rate = fitg$estimate[2]),
  # #                            color = Colors.VIOLET);
  # #
  # # logNormalPlotAndHist <- Plotter.getLogNormalDistributionPlot(scores, fitln, Colors.BLUE, print = FALSE);
  # #
  # # plot <- logNormalPlotAndHist + weibullPlot + gammaPlot;
  # #
  # # print(plot);
  # 
  # data <- fitdistrplus::gofstat(list(fitW, fitg, fitln), fitnames = c("weibull", "gamma", "lognormal"))
  # print(as.character(data$kstest["weibull"]) == "not rejected");
  
  cumulativeDistributionFunction <- plnorm;

  # find out correct distribution parameters
  fit <- fitdistrplus::fitdist(scores, "lnorm");  # lognormal distribution

  # compute p-value under lognormal distributions
  sortedScores <- sort(scores);
  pValue <- cumulativeDistributionFunction(sortedScores[[threshold]], fit$estimate[1], fit$estimate[2]);  # P[X <= threshold]

  # hist(scores, pch=20, breaks=25, prob=TRUE, main = Symbols.EMPTY);
  # curve(densityFunction(x, fit$estimate[1], fit$estimate[2]), col = "blue", lwd = 2, add = T);

  # draw
  # Plotter.getLogNormalDistributionPlot(scores, fit, Colors.BLUE, print = TRUE);

  return(pValue);
}


#########################################################
#' Returns ranks and corresponding p-vlaues.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreTypes {list} contains different scoreTypes
#' @param threshold {numeric} the threshold for which the p-value has to be computed
#' @param computeDeltaScores {logical} tells if the p-value for delta-scores should be computed
#'
#' @return {list} p-values and corresponding ranks
Analyzer.getPValuesAndCorrespondingRanks <- function(passPathGenericName, scoresFolderGenericName, 
                                                     passes, scoreTypes, threshold, computeDeltaScores = FALSE) {
  
  # compute corrects scores ranks and p-values of the scores
  pValuesPerPass <- list();
  ranksPerPass <- list();
  
  # iterate over the different passes
  for (i in 1:length(passes)) {
    passesPath <- paste(passPathGenericName, passes[[i]], sep = Symbols.EMPTY);
    
    ranksPerScoreType <- list();
    pValuesPerScoreType <- list();
    
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
      
      ranks <- Analyzer.getRanks(scoresForCorrectPositions, scoresPerSample, FALSE);
      pValues <- Analyzer.getPValues(scoresPerSample, threshold, computeDeltaScores);
      
      # store the ranks per score-type
      ranksPerScoreType <- rlist::list.append(ranksPerScoreType, ranks);
      pValuesPerScoreType <- rlist::list.append(pValuesPerScoreType, pValues);
    }
    
    ranksPerPass <- rlist::list.append(ranksPerPass, ranksPerScoreType);
    pValuesPerPass <- rlist::list.append(pValuesPerPass, pValuesPerScoreType)
  }
  
  return(list(ranksPerPass = ranksPerPass, pValuesPerPass = pValuesPerPass));
}


#########################################################
#' Returns the profile index with the lowest pairwise distance.
#'
#' @param cluster {list} the cluster from which the profile with the lowest pairwise distance should be retuned
#'
#' @return {list} the profile with lowest pairwise ditance
#' @export
Analyzer.__getIndexOfMostSimilarProfile <- function(cluster) {
  if (length(cluster) == 1) return(1);
  
  distances <- list();
  
  # iterate over each profile in the cluster
  for (i in 1:length(cluster)) {
    profile <- cluster[[i]];
    remainingProfiles <- cluster[-i];
    
    pairwiseDistancesSum <- 0;
    
    for (j in seq_len(length(remainingProfiles))) {
      secondProfile <- remainingProfiles[[j]];
      numberOfPoints <- length(profile);
      distance <- sum(abs(profile - secondProfile))/numberOfPoints;
      pairwiseDistancesSum <- pairwiseDistancesSum + distance;
    }
    
    distances <- rlist::list.append(distances, pairwiseDistancesSum);
  }
  
  # find profile with minimum distance
  minDistance <- min(unlist(distances));
  bestPosition <- which(distances == minDistance)[[1]];
  
  return(bestPosition);
}


#########################################################
#' Sorts the scores in ascending order and computes differences to the best score. 
#' These dfferences are returned.
#'
#' @param scoresOfSample {list} the scores created by shifting a single sample in a chronology and computing distances
#'
#' @return {list} the difference scores
#' @export
Analyzer.getDeltaScores <- function(scoresOfSample) {
  scores <- unlist(scoresOfSample);
  
  # sort scores
  sortedScores <- sort(scores);
  
  # find best score
  bestScore <- sortedScores[[1]];
  remainingScores <- sortedScores[-1];
  
  return(remainingScores-bestScore);
}


#########################################################
#' Computes for each sample the counts for each position.
#'
#' @param bucketChronology {list} the list of buckets
#' @param namesChronology {list} the list of per bucket used trees
#' @param samples {list} data of multiple curve-files (profiles with its years)
#'
#' @return {list} the counts per year for each sample
#' @export
Analyzer.computeCountsPerSample <- function(bucketChronology, namesChronology, samples) {
  countsPerSample <- list();
  
  # prepare data
  data <- Analyzer.__dataPreparation2(bucketChronology, samples, TRUE, TRUE);
  namesChronology <- Analyzer.__EnvToListChronology(namesChronology);
  
  print(Titles.CURRENT_PROCESSED_SAMPLE);
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    patternProfiles <- data$samples[[i]]$profiles;
    
    counts <- Analyzer.__computeTreeCounts(data$bucketChronology, namesChronology, patternProfiles);
    
    countsPerSample <- rlist::list.append(countsPerSample, counts);
  }
  
  return(countsPerSample);
}


#########################################################
#' Reconverts an environment chronology into a list-based chronology.
#'
#' @param chronology {list} the environment which should be converted
#'
#' @return {list} list-based chronology
Analyzer.__EnvToListChronology <- function(chronology) {
  editedChronology <- list();
  
  # iterate over all years and compute for every profile in a year the desired operations
  for (i in ls(chronology)) {
    bucket <- chronology[[i]];
    
    if (!is.null(bucket)) {  # if bucket exist
      editedBucket <- list();
      
      # iterate over all data in bucket
      for (j in 1:length(bucket)) {
        data <- bucket[[j]];
        
        editedBucket <- rlist::list.append(editedBucket, data);
      }
      
      editedChronology <- rlist::list.append(editedChronology, editedBucket);
    }
  }
  
  return(editedChronology);
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the bucket-chronology, which consists out of buckets of profiles
#' and computes the number of trees used for the final distance.
#'
#' @param consensusProfiles {list} the list of profiles of the consensus
#' @param namesChronology {list} the list of per bucket used trees
#' @param patternProfiles {list} the list of profiles of the pattern
#' 
#' @return {list} the scores
Analyzer.__computeTreeCounts <- function(bucketChronology, namesChronology, patternProfiles) {
  counts <- list();
  
  n <- length(bucketChronology);
  m <- length(patternProfiles);
  
  # iterate over each curve position with the subpattern
  for (i in 1:(n-m+1)) {  # "+1", because i is "1" at the beginning and else we wouldn't align against every position
    # grab a cuts of pattern-size
    bucketCut <- CurvesMiner.getSubset(bucketChronology, i, (i-1) + m);  # "-1", because i is "1" at the beginning
    namesCut <- CurvesMiner.getSubset(namesChronology, i, (i-1) + m);
    
    accumulatedScore <- 0;
    namesPerSampleProfile <- list();
    
    # iterate over each position of the pattern and the cuts
    for (k in 1:m) {
      bucket <- bucketCut[[k]];
      names <- namesCut[[k]];
      subpatternProfile <- patternProfiles[[k]];
      
      bucketScores <- list();
      treeNames <- list();
      
      # iterate over all profiles and tree names in the bucket
      for (l in 1:length(bucket)) {
        bucketProfile <- bucket[[l]];
        profileTreeName <- names[[l]];
        
        alignmentScore <- Analyzer.__getScore(subpatternProfile, bucketProfile);
        
        bucketScores <- rlist::list.append(bucketScores, alignmentScore);
        treeNames <- rlist::list.append(treeNames, profileTreeName);
      }
      
      bestScore <- min(unlist(bucketScores));
      bestScoreIndex <- which(bucketScores == bestScore);
      bestScoreName <- treeNames[[bestScoreIndex]];
      
      namesPerSampleProfile <- rlist::list.append(namesPerSampleProfile, bestScoreName);
    }
    
    count <- unique(unlist(namesPerSampleProfile));
    
    counts <- rlist::list.append(counts, length(count)); 
  }
  
  return(counts);
} 


#########################################################
#' Extracts the score of each sample at its correct position in the consensus.
#'
#' @param countsPerSample {list} the scores per sample
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param years {list} the years of the consensus
#'
#' @return {list} score per sample
Analyzer.getCountsForCorrectPositions <- function(countsPerSample, samples, years) {
  counts <- list();
  
  # get for each sample the correct score
  for (i in 1:length(samples)) {
    countsOfSample <- countsPerSample[[i]];
    sample <- samples[[i]];

    correctCount <- Analyzer.__getCorrectPatternPositionData(countsOfSample, 
                                                             sample$years[[1]], years);
    counts <- rlist::list.append(counts, correctCount);
  }
  
  return(counts);
}


#########################################################
#' Returns the start- and end-years out of the per tree data structure.
#'
#' @param perTreeDataStructure {list} the per tree data-structure
#'
#' @return {list(startYears, endYears)} the list containing lists of start- and end-years
#' @export
Analyzer.getStartAndEndYears <- function(perTreeDataStructure) {
  startYears <- list();
  endYears <- list();
  
  for (i in seq_along(perTreeDataStructure)) {
    tree <- perTreeDataStructure[[i]];
    
    startYear <- tree$years[[1]];
    endYear <- tree$years[[length(tree$years)]];
    
    startYears <- rlist::list.append(startYears, startYear);
    endYears <- rlist::list.append(endYears, endYear);
  }
  
  return(list(startYears = startYears, endYears = endYears));
}


#########################################################
#' Computes the scores for the per-tree approach.
#'
#' @param perTreeChronology {list} the list of trees
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#' @param startYears {list} the years at which the chronologies starts
#' @param endYears {list} the years at which the chronologies ends
#' @param yearsToLookAt {list} the chronology years
#' 
#' @return {list} the scores per sample
#' @export
Analyzer.computeScoresPerSample4 <- function(perTreeChronology, samples, normalize, slopeBased,
                                             startYears, endYears, yearsToLookAt) {
  scoresPerSample <- list();
  
  # prepare data
  data <- Analyzer.__dataPreparation3(perTreeChronology, samples, normalize, slopeBased);
  perTreeChronology <- data$perTreeChronology;
  
  print(Titles.CURRENT_PROCESSED_SAMPLE);
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    patternProfiles <- data$samples[[i]]$profiles;
    
    scores <- Analyzer.__computePerTreeScores(perTreeChronology, patternProfiles, startYears, endYears, yearsToLookAt);
    
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Prepares the data for the per tree approach i.e. replaces with the slopes or normalized data.
#'
#' @param perTreeChronology {envrionment} a chronology that stores all sequences seperatly
#' @param samples {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param slopeBased {logical} this parameter can be set to compute score on slopes
#'
#' @return {list(perTreeChronology, samples)} the prepared data
Analyzer.__dataPreparation3 <- function(perTreeChronology, samples, normalize, slopeBased) {
  
  editedProfiles <- Math.computeOperation(perTreeChronology, normalize, slopeBased);
  
  # iterate over each tree
  for (i in seq_along(perTreeChronology)) {
    tree <- perTreeChronology[[i]];
    profiles <- tree$profiles;

    # write back
    tree$profiles <- editedProfiles[[i]][[1]];
    perTreeChronology[[i]] <- tree;
  }

  samples <- Math.computeOperation(samples, normalize, slopeBased);
  
  return(list(perTreeChronology = perTreeChronology, samples = samples));
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the bucket-chronology, which consists out of buckets of profiles.
#'
#' @param perTreeChronology {list} the list of profiles of the consensus
#' @param patternProfiles {list} the list of profiles of the pattern
#' @param startYears {list} the years at which the chronologies starts
#' @param endYears {list} the years at which the chronologies ends
#' @param yearsToLookAt {list} the chronology years
#' 
#' @return {list} the scores
Analyzer.__computePerTreeScores <- function(perTreeChronology, patternProfiles, 
                                            startYears, endYears, yearsToLookAt) {
  scores <- list();
  
  n <- length(yearsToLookAt);
  m <- length(patternProfiles);
  
  # init positions in chronologies with ones
  currentPositionsInChronologies <- rep(1, length(perTreeChronology));  # stores the current local positions per chronology

  # iterate over all years at which you have to look at
  for (i in 1:(n-m+1)) {
    year <- yearsToLookAt[[i]];
    
    activeStartYearIndices <- which(startYears <= year);  # bigger/equal left limit
    activeEndYearIndices <- which(endYears >= year);  # lower/equal right limit
    activeTreeIndices <- intersect(activeStartYearIndices, activeEndYearIndices);
    
    # extract active trees
    activeTrees <- perTreeChronology[activeTreeIndices];
    
    perTreeScores <- list();
    
    # iterate over active trees
    for (j in seq_along(activeTrees)) {
      # extract active tree
      tree <- activeTrees[[j]];
      treeYears <- tree$years;
      treeProfiles <- tree$profiles;
      
      # extract current position in chronology
      activeIndex <- activeTreeIndices[[j]];
      currentChronologyPosition <- currentPositionsInChronologies[[activeIndex]];
      
      startYear <- as.numeric(treeYears[[currentChronologyPosition]]);
      
      # check if we are not already further in the chronology
      if (startYear == year) {
        
        # extract from tree if possible
        if (currentChronologyPosition + m - 1 <= length(treeYears)) {
  
          endYear <- as.numeric(treeYears[[currentChronologyPosition + m - 1]]);
  
            if (endYear == year + m - 1) {  # if consecutive years without gaps
              # extract profiles for corresponding years
              cut <- CurvesMiner.getSubset(treeProfiles, currentChronologyPosition, (currentChronologyPosition-1) + m);
              
              sampleScore <- 0;
              
              # iterate over each position of the sample and the cut
              for (k in 1:length(cut)) {
                profile <- cut[[k]];
                sampleProfile <- patternProfiles[[k]];
                
                if (identical(profile, sampleProfile)) { 
                  print(Strings.BUG);  # this should be impossible
                  View(profile);
                  View(sampleProfile);
                  stop("DONE");
                }
                
                alignmentScore <- Analyzer.__getScore(profile, sampleProfile);
                
                sampleScore <- sampleScore + alignmentScore;
              }
              
              perTreeScores <- rlist::list.append(perTreeScores, sampleScore);
            }
        } 
        
        # increment "activeIndex" location in "currentPositionsInChronologies"-vector (incrementing counter for current chronology)
          currentPositionsInChronologies[[activeIndex]] <- currentChronologyPosition + 1;
      }
    }
    
    aggregatedScores <- unlist(perTreeScores);
    
    if (length(aggregatedScores) == 0) {
      scores <- rlist::list.append(scores, Inf);
    } else {
      scores <- rlist::list.append(scores, min(aggregatedScores));  # selects from all tree-scores, the lowest
    }
  }
  
  return(scores);
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the bucket-chronology, which consists out of buckets of profiles.
#'
#' @param chronology {list} the list of buckets
#' @param sample {list} the list of profiles of the pattern
#' @param perYearScores {list} the list in which it is stored
#' 
#' @return {list} the per year scores
Analyzer.computeScoresOfSample <- function(chronology, sample, perYearScores) {
  years <- chronology$years;
  
  n <- length(years);
  m <- length(sample);

  # iterate over all years
  if (n >= m) {
    for (i in 1:(n-m+1)) {
      startYear <- as.numeric(years[[i]]);
      endYear <- as.numeric(years[[i + m - 1]]);
      
      # extract profiles for corresponding years
      cut <- CurvesMiner.getSubset(chronology$profiles, i, (i-1) + m);
      
      # check if consecutive sample can be cut out
      if (endYear - startYear == m - 1) {
        # compute score
        
        score <- 0;
        
        for (k in 1:length(cut)) {
          cutProfile <- unlist(cut[[k]]);
          profile <- unlist(sample[[k]]);
          
          alignmentScore <- Analyzer.__getScore(profile, cutProfile);
          
          score <- score + alignmentScore;
        }
        
        # store score in right year
        if (is.null(perYearScores[[toString(startYear)]])) {
          perYearScores[[toString(startYear)]] <- list(score);
        } else {
          perYearScores[[toString(startYear)]] <- rlist::list.append(perYearScores[[toString(startYear)]], score);
        }
      }
    }
  }
  
  return(perYearScores);
}


#########################################################
#' Computes correlation coefficients 
#' at each position of the sample within each chronology.
#'
#' @param chronologies {list} the list of per-tree chronologies
#' @param samples {list} data of multiple curve-files (widths with its years)
#' @param method {string} the correlation method which should be used {"spearman", "pearson", "kendall"}
#' @param startYears {list} the years at which the chronologies starts
#' @param endYears {list} the years at which the chronologies ends
#' @param yearsToLookAt {list} the chronology years
#' 
#' @return {list} the scores per sample
#' @export
Analyzer.computeCoefficientsPerSample2 <- function(chronologies, samples, method, startYears, endYears, yearsToLookAt) {
  scoresPerSample <- list();
  
  # iterate over all samples
  for (i in 1:length(samples)) {
    print(i);
    widths <- samples[[i]];
    scores <- Analyzer.__computeCoefficients2(chronologies, unlist(widths), method,
                                              startYears, endYears, yearsToLookAt);
    scoresPerSample <- rlist::list.append(scoresPerSample, scores);
  }
  
  return(scoresPerSample);
}


#########################################################
#' Aligns every pattern consisting out of profiles against every position
#' in the consensus, which consists out of profiles.
#'
#' @param perTreeChronology {list} the list of trees
#' @param pattern {vector} the pattern which have to be found in the chronology
#' @param method {string} the correlation coefficient method
#' @param startYears {list} the years at which the chronologies starts
#' @param endYears {list} the years at which the chronologies ends
#' @param yearsToLookAt {list} the chronology years
#'
#' @return {list} the scores per sample
#' @export
Analyzer.__computeCoefficients2 <- function(perTreeChronology, pattern, method, 
                                            startYears, endYears, yearsToLookAt) {
  scores <- list();
  
  n <- length(yearsToLookAt);
  m <- length(pattern);
  
  # init positions in chronologies with ones
  currentPositionsInChronologies <- rep(1, length(perTreeChronology));  # stores the current local positions per chronology
  
  # iterate over all years at which you have to look at
  for (i in 1:(n-m+1)) {
    year <- yearsToLookAt[[i]];
    
    activeStartYearIndices <- which(startYears <= year);  # bigger/equal left limit
    activeEndYearIndices <- which(endYears >= year);  # lower/equal right limit
    activeTreeIndices <- intersect(activeStartYearIndices, activeEndYearIndices);
    
    # extract active trees
    activeTrees <- perTreeChronology[activeTreeIndices];
    
    perTreeScores <- list();
    
    # iterate over active trees
    for (j in seq_along(activeTrees)) {
      # extract active tree
      tree <- activeTrees[[j]];
      treeYears <- tree$years;
      treeWidths <- tree$widths;
      
      # extract current position in chronology
      activeIndex <- activeTreeIndices[[j]];
      currentChronologyPosition <- currentPositionsInChronologies[[activeIndex]];
      
      startYear <- as.numeric(treeYears[[currentChronologyPosition]]);
      
      # check if we are not already further in the chronology
      if (startYear == year) {
        
        # extract from tree if possible
        if (currentChronologyPosition + m - 1 <= length(treeYears)) {
          
          endYear <- as.numeric(treeYears[[currentChronologyPosition + m - 1]]);
          
          if (endYear == year + m - 1) {  # if consecutive years without gaps
            # extract widths for corresponding years
            cut <- CurvesMiner.getSubset(treeWidths, currentChronologyPosition, (currentChronologyPosition-1) + m);
            
            if (identical(method, Defaults.CORRELATION_T_VALUE)) {
              sampleScore <- Math.twoSampleTtestValues(unlist(cut), pattern);
            } else {
              sampleScore <- cor(unlist(cut), pattern, method = method);
            }
            
            perTreeScores <- rlist::list.append(perTreeScores, sampleScore);
          }
        } 
        
        # increment "activeIndex" location in "currentPositionsInChronologies"-vector (incrementing counter for current chronology)
        currentPositionsInChronologies[[activeIndex]] <- currentChronologyPosition + 1;
      }
    }
    
    aggregatedScores <- unlist(perTreeScores);
    
    if (length(aggregatedScores) == 0) {
      scores <- rlist::list.append(scores, -Inf);
    } else {
      scores <- rlist::list.append(scores, max(aggregatedScores));  # selects from all tree-scores, the highest
    }
  }
  
  return(scores);
}
