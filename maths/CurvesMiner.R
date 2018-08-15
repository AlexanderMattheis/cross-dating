#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

CURVE_MINER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Modifies, filters and adds/concatenates curve data.

#########################################################
#' Filters the data of a single tree.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#' @param treeNumber {numerical} the tree from which you want get back data
#'
#' @return {list} curves of a single tree in all measured directions
#' @export
CurvesMiner.filterTree <- function(curvesData, treeNumber) {
  numberOfCurves <- length(curvesData);
  
  curvesWithGivenTreeNumber <- list();
  
  # iterate over all curves and filter out all with the given treeNumber
  for (i in 1:numberOfCurves) {
    curveData <- curvesData[[i]];
    curveName <- curveData$name;
    curveTreeNumber <- as.numeric(substr(stringr::str_extract_all(curveName, Expression.NUMBERS), 1, 2)); 
    
    if (curveTreeNumber == treeNumber) {
      curvesWithGivenTreeNumber <- rlist::list.append(curvesWithGivenTreeNumber, curveData);
    }
  }
  
  return(curvesWithGivenTreeNumber);
}


#########################################################
#' Computes profile-wise consensus of all directions.
#' Hint: The curves must have exactly the same years without wrong values.
#'
#' @param treeData {list} data from several curve-files of a single tree (profiles with its years)
#' @param distanceFunction {numerical} the distance function
#' 
#' @return {list} curves of a single tree in all measured directions
#' @export
CurvesMiner.computeTreeConsensus <- function(treeData, distanceFunction) {
  numberOfDirections <- length(treeData);
  consensusProfilesList <- list();
  
  if (numberOfDirections > 0) {
    firstCurve <- treeData[[1]];
    numberOfYears <- length(firstCurve$years);  # hint: all curves should have same years and number of years
    
    # iterate over all profiles (for each year)
    for (i in 1:numberOfYears) {
      profilesToAlignList <- list();
      
      # iterate over all diretions
      for (j in 1:numberOfDirections) {
        profilesToAlignList <- rlist::list.append(profilesToAlignList, treeData[[j]]$profiles[[i]]);
      }
      
      profilesToAlign <- CurvesMiner.__equalizeAllPointNumbers(profilesToAlignList);
      consensusProfileData <- alignCurves(y = profilesToAlign, distFunc = distanceFunction);
      consensusProfilesList <- rlist::list.append(consensusProfilesList, consensusProfileData$consensus$y);
    }
  }
  
  return(consensusProfilesList);
}


#########################################################
#' Adapts profiles to an equal number of points
#' by adding Not-Availables at the end of too short profiles.
#'
#' @param profiles {list} multiple profiles which should be equalized
#'
#' @return {dataframe} equalized profiles
CurvesMiner.__equalizeAllPointNumbers <- function(profiles) {
  maxLength <- -1;
  
  # iterate over all profiles to find out the max number of points in a profile
  for (i in 1:length(profiles)) {
    profile <- profiles[[i]];
    numOfProfilePoints <- length(profile);
    
    if (maxLength < numOfProfilePoints) {
      maxLength <- numOfProfilePoints;
    }
  }
  
  # iterate over all profiles to equalize each profile to same number of points
  equalizedProfiles <- list();
  for (i in 1:length(profiles)) {
    profile <- unlist(profiles[[i]]);
    equalizedProfile <- c(profile, rep(NA, maxLength - length(profile)));
    
    equalizedProfiles <- rlist::list.append(equalizedProfiles, equalizedProfile);
  }
  
  dataframe <- (do.call(cbind.data.frame, equalizedProfiles));
  names(dataframe) <- 1:ncol(dataframe);  # names are necessary
  return(dataframe);
}


#########################################################
#' Subsets a profiles list.
#'
#' @param profiles {list} the set from which you want a subset
#' @param start {numerical} the first position from the list which should be used
#' @param end {numerical} the last position from thelist which should be used
#'
#' @return {list} a subset
#' @export
CurvesMiner.getSubset <- function(set, start, end) {
  subset <- list();
  
  if (start >= 1 && end <= length(set)) {
    for (i in start:end) {
      subset <- rlist::list.append(subset, set[[i]]);
    }
  }

  return(subset);
}


#########################################################
#' Extracts unifomly a number of samples from curve data.
#' Hint: It is possible that there is not enough data or 
#' that two times the same years and file is chosen. 
#' In this case you get not your desired number of samples.
#'
#' @param numberOfSamples {numerical} the desired number of samples
#' @param sampleLength {numerical} the lenght of the samples which are extracted
#' @param curvesData {list} data of multiple curve-files (profiles with their years)
#' @param widthBased {logical} tells if widths should be used in instead of profiles
#'
#' @return {list(samples, years, names)} list of extracted samples
#' @export
CurvesMiner.extractUniformlyConsecutiveSamples <- function(numberOfSamples, sampleLength, curvesData, widthBased = FALSE) {
  sampleList <- CurvesMiner.__createSampleList(curvesData, sampleLength, widthBased);
  return(CurvesMiner.__selectSamples(sampleList, numberOfSamples, sampleLength));
}


#########################################################
#' Creates a list of samples with conscecutive years i.e. samples with 
#' non-consecutive years are sorted out.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with their years)
#' @param sampleLength {numerical} the lenght of the samples which are extracted
#' @param widthBased {logical} tells if widths should be used in instead of profiles
#' 
#' @return {list(years, samples, names)} the data to select from
#' @example 
#' year | profile   | name
#' 1992   p1 p2 p3   "401_MICA-cons"
#' 1993   p2 p3 p4   "401_MICA-cons"
CurvesMiner.__createSampleList <- function(curvesData, sampleLength, widthBased) {
  yearsList <- list();
  samplesList <- list();
  treeNamesList <- list();
  
  # iterate over all curves
  for (i in 1:length(curvesData)) {  # iterate over all curves
    # read in data
    curve <- curvesData[[i]];
    
    if (length(curve$years) >= sampleLength) {
      # iterate over each year
      for (j in 1:(length(curve$years) - sampleLength + 1)) {
        # create sample of length sampleLength
        sampleStartYear <- curve$years[[j]];
        sampleEndYear <- curve$years[[j + sampleLength - 1]];
        
        if (as.numeric(sampleStartYear) + sampleLength - 1 == as.numeric(sampleEndYear)) {  # only if consecutive years
          if (widthBased) {
            sample <- curve$widths[j:(j + sampleLength - 1)];
          } else {
            sample <- curve$profiles[j:(j + sampleLength - 1)];
          }

          treeName <- curve$name;
          
          yearsList <- rlist::list.append(yearsList, sampleStartYear);
          samplesList <- rlist::list.append(samplesList, sample);
          treeNamesList <- rlist::list.append(treeNamesList, treeName);
        }
      }
    }
  }
   
  return(list(years = yearsList, samples = samplesList, names = treeNamesList));
}


#########################################################
#' Selects uniformly samples.
#' Hint: It tries to select samples. But if there
#' are not enough, then it can happen, that not the desired number of samples
#' is returned, but a lower one.
#'
#' @param sampleList {list(years, profiles|widths, names)} converted data of multiple curve-files 
#' (every profile with its year and curvename)
#' @param numberOfSamples {numerical} the desired number of samples
#' @param sampleLength {numerical} the length of the samples which are extracteds
#' 
#' @return {list(samples, years, names)} the samples with its years and names seperatly
CurvesMiner.__selectSamples <- function(sampleList, numberOfSamples, sampleLength) {
  # to create
  samples <- list();
  samplesYears <- list();
  curveNames <- list();
  
  # to avoid selecting two times the same index
  selectedIndices <- list();  
  
  # extract data
  points <- sampleList$samples;
  years <- sampleList$years;
  names <- sampleList$names;

  # create samples
  numOfYearIndices <- length(sampleList$years);  # 1992 1993, 1990 .. from multiple files
  print(Strings.PROFILES);
  print(numOfYearIndices);
  tries <- 0;
  
  while (length(samples) < numberOfSamples && 
         tries < Parameters.TRIES_TO_FIND_A_SAMPLE) { 
    indexToSelect <- sample(1:numOfYearIndices, 1, replace = TRUE);
    
    tries <- tries + 1;
    
    if (!(indexToSelect %in% selectedIndices)) {
      samples <- rlist::list.append(samples, points[[indexToSelect]]);
      samplesYears <- rlist::list.append(samplesYears, years[[indexToSelect]]);
      curveNames <- rlist::list.append(curveNames, names[[indexToSelect]]);
      
      selectedIndices <- rlist::list.append(selectedIndices, indexToSelect);
      tries <- 0;
    }
  }
  
  return(list(samples = samples, years = samplesYears, names = curveNames));
}


#########################################################
#' Removes all rows from the sample data in the curve data.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#' @param sampleData {list(samples, years, curveNames)} list of extracted samples with consecutive years
#' @param widthBased {logical} tells if widths should be used in instead of profiles
#' 
#' @return {list} curveData in which the rows from the sampleData has been removed
#' @export
CurvesMiner.subtractSampleData <- function(curvesData, sampleData, widthBased = FALSE) {
  indicesToRemove <- CurvesMiner.__computeIndicesToRemove(curvesData, sampleData);
  
  for (i in 1:length(curvesData)) {
    if(!is.null(indicesToRemove[[toString(i)]])) {
      toRemove <- unique(indicesToRemove[[toString(i)]]);  
      toRemove <- toRemove[1:length(toRemove)];  # avoids a NULL at the end
      
      if (widthBased) {  # extension for ring widths
        curvesData[[i]]$widths <- curvesData[[i]]$widths[-toRemove];
      } else {
        curvesData[[i]]$profiles <- curvesData[[i]]$profiles[-toRemove];
      }

      curvesData[[i]]$years <- curvesData[[i]]$years[-toRemove];
    }
  }
  
  # View(sampleData);
  # View(curvesData);
  
  return(curvesData);
}


#########################################################
#' Computes the indices to remove from curvesData.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#' @param sampleData {list(samples, years, curveNames)} list of extracted samples with consecutive years
#'
#' @return {environment("numbers", vectors)} the indices which should be removed from curvesData
CurvesMiner.__computeIndicesToRemove <- function(curvesData, sampleData) {
  toRemoveIndices <- new.env();
  
  # iterate over all samples in sampleData
  for (i in 1:length(sampleData$samples)) {
    # look up name of curve from which the sample has been extracted
    curveFileName <- sampleData$names[[i]];
    
    # look up the first year, which was removed
    curveYear <- sampleData$years[[i]];
    
    # look up how many consecutive years should be removed: 
    # works only for consecutive years (but this is given with the CurvesMiner.extractUniformlyConsecutiveSamples function)
    numOfYears <- length(sampleData$samples[[i]]);
    
    # iterate over all curve files, to find that curveFileName
    for (j in 1:length(curvesData)) {
      curve <- curvesData[[j]];
      
      if (curve$name == curveFileName) {  # if curve has been found
        
        # iterate over all years/profiles in that file
        for (k in 1:length(curve$years)) {
          year <- curve$years[[k]];
          
          if (year == curveYear) {  # if startyear for removement has been found
            # add year/profile indices which should be removed to a list
            
            if(is.null(toRemoveIndices[[toString(j)]])) {  # if filled for the first time
              toRemoveIndices[[toString(j)]] <- k:(k + numOfYears - 1);  # "-1", because k is one at begin.
            } else {  # if has been filled before
              toRemoveIndices[[toString(j)]] <- c(toRemoveIndices[[toString(j)]], k:(k + numOfYears - 1));
            }
            
            # print(i);  # test
            found <- TRUE;  # mark, to jump out of the mid loop
            break;  # jump out of the most inner loop
          }
          
        }
        
        if (found == TRUE) break;
      }
    }
  }

  return(toRemoveIndices);
}


#########################################################
#' Collect per year all widths to compute later a consensus for these years.
#' It is analogous to the CurvesMiner.collectData function for profiles.
#' But without less checkups since the widths will be also compared 
#' with the profiles which are checked.
#' Sadly the widths were received two months after the profiles such that
#' a generic function couldn't be written anymore since the dependencies
#' in the project were no more reconstructable due to the huge amount of code.
#' 
#' @param curvesData {list} data of multiple curve-files (points with its years)
#' 
#' @return {list(environment("years", "list of profiles"), environment("years", "list of names"), minStartYear, maxEndYear)} 
#' the profiles & names per year, as well as the minimum start- and maximum end-year
#' @export
CurvesMiner.collectWidthsData <- function(curvesData) {
  widthsPerYear <- new.env();
  namesPerYear <- new.env();
  minStartYear <- .Machine$integer.max;
  maxEndYear <- -1;
  
  # iterate over all curves
  for (i in 1:length(curvesData)) {
    curve <- curvesData[[i]];
    years <- curve$years;
    widths <- curve$widths;
    
    if (is.null(widths)) {
      widths <- curve$characteristics;
    }
    
    name <- curve$name;
    
    # iterate over all years
    for (j in 1:length(years)) {
      year <- years[[j]];
      width <- widths[[j]];  # assigned profile
      
      if (is.null(widthsPerYear[[toString(year)]])) {
        widthsPerYear[[toString(year)]] <- list(width);
        namesPerYear[[toString(year)]] <- list(name);
      } else {
        widthsPerYear[[toString(year)]] <- rlist::list.append(widthsPerYear[[toString(year)]], width);
        namesPerYear[[toString(year)]] <- rlist::list.append(namesPerYear[[toString(year)]], name);
      }
    }
    
    if (length(years) > 0) {  # just for the case the user inputs something what is far too short
      if (minStartYear > years[[1]]) {
        minStartYear <- years[[1]];
      }
      
      if (maxEndYear < years[[length(years)]]) {
        maxEndYear <- years[[length(years)]];
      }
    }
  }
  
  return(list(widthsPerYear = widthsPerYear, namesPerYear = namesPerYear, 
              minStartYear = minStartYear, maxEndYear = maxEndYear));
}


#########################################################
#' Collect per year all profiles to compute later a consensus for these years.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#'
#' @return {environment("years", "list of profiles"), environment("years", "list of names")} the profiles per year
#' @export
CurvesMiner.collectData <- function(curvesData) {
  profilesPerYear <- new.env();
  namesPerYear <- new.env();
  minStartYear <- .Machine$integer.max;
  maxStartYear <- -1;
  
  # iterate over all curves
  for (i in 1:length(curvesData)) {
    curve <- curvesData[[i]];
    years <- curve$years;
    profiles <- curve$profiles;
    name <- curve$name;
    
    # iterate over all years
    for (j in 1:length(years)) {
      year <- years[[j]];
      profile <- profiles[[j]];  # assigned profile
      
      if (is.null(profilesPerYear[[toString(year)]])) {
        profilesPerYear[[toString(year)]] <- list(profile);
        namesPerYear[[toString(year)]] <- list(name);
      } else {
        if (!CurvesMiner.__contained(profile, profilesPerYear[[toString(year)]], name)) {
          profilesPerYear[[toString(year)]] <- rlist::list.append(profilesPerYear[[toString(year)]], profile);
          namesPerYear[[toString(year)]] <- rlist::list.append(namesPerYear[[toString(year)]], name);
        }
      }
    }
    
    if (length(years) > 0) {  # just for the case the user inputs something what is far too short
      if (minStartYear > years[[1]]) {
        minStartYear <- years[[1]];
      }
      
      if (maxStartYear < years[[length(years)]]) {
        maxStartYear <- years[[length(years)]];
      }
    }
  }

  numberOfProfiles <- 0;
  
  # check for gaps
  print(Titles.NUMBER_OF_PROFILES);
  for (i in ls(profilesPerYear)) {
    print(i);
    numberOfProfiles <- numberOfProfiles + length(profilesPerYear[[i]]);
    print(length(profilesPerYear[[i]]));
    print(Strings.SEPERATION);
  }
  
  print(Titles.NUMBER_OF_PROFILES_FINAL);
  print(numberOfProfiles);

  return(list(profilesPerYear = profilesPerYear, namesPerYear = namesPerYear, 
              minStartYear = minStartYear, maxStartYear = maxStartYear));
}


#########################################################
#' Checks if an element is contained in another list.
#'
#' @param elementToFind {list} the element which has to be find
#' @param set {list} the set in which the element has to be searched
#' @param name {string} the name of the file
#'
#' @return {logical} TRUE or FALSE depending on if the element is in the set or not
CurvesMiner.__contained <- function(elementToFind, set, name) {
  for (i in 1:length(set)) {
    element <- set[[i]];
    
    if (identical(element, elementToFind)) {
      print(name);
      print(Strings.DATASET_CORRUPT);
      return(TRUE);
    }
  }

  
  return(FALSE);
}


#########################################################
#' Computes the consensus per year.
#'
#' @param profilesPerYear {environment("years", "list of profiles")} the profiles per year
#' @param startYear {numerical} the year in which the master-chronology should start
#' @param endYear {numerical} the year in which the master-chronology should end
#' @param plot {logical} tells if a plot should be created or not
#' @param yLimits {vector} the limits for the plot
#'
#' @return {list(years, profiles, numberOfProfiles)} the unique years and the profiles
#' @export
CurvesMiner.computeChronology <- function(profilesPerYear, startYear, endYear, plot, yLimits) {
  consensusProfilesList <- list();
  years <- list();
  numberOfProfilesList <- list();
  
  addedToProfileNumberList <- FALSE;
  
  print(Titles.START_CONSENSUS_COMPUTATION);
  
  # iterate over all years
  for (i in startYear:endYear) {
    addedToProfileNumberList <- FALSE;
    
    print(i); # to see if it is still working
    
    profiles <- profilesPerYear[[toString(i)]];
    
    if (!is.null(profiles)) {  # if for the year profiles exist
      # subtract mean and divide by standard deviation
      normalizedProfiles <- Math.computePureProfiles(profiles);

      numProfiles <- length(normalizedProfiles);
      
      # fixing bug in MICA
      if (numProfiles >= Parameters.LIMIT_PROFILES_PER_CONSENSUS) {  # delete some columns or the alignment will fail or take too much time
        data <- Clustering.getSubsetOfProfiles2(normalizedProfiles, profiles);
        profilesToAlign <- CurvesMiner.__equalizeAllPointNumbers(data$normalizedProfiles);  # the profiles which are aligned
        profilesForAlignment <- CurvesMiner.__equalizeAllPointNumbers(data$profiles);  # the profiles on which alignment is applied 
        numProfiles <- ncol(profilesToAlign);
      } else {
        profilesToAlign <- CurvesMiner.__equalizeAllPointNumbers(normalizedProfiles);  # the profiles which are aligned
        profilesForAlignment <- CurvesMiner.__equalizeAllPointNumbers(profiles);  # the profiles on which alignment is applied 
      }
      
      if (numProfiles > 1) {  # if there is something to align
        consensusProfileData <- alignCurves(y = profilesToAlign, 
                                            distFunc = Parameters.DIST_FUNCTION,
                                            distSample = Parameters.DIST_SAMPLE,
                                            distWarpScaling = Parameters.DIST_WARP_SCALING,
                                            maxWarpingFactor = Parameters.MAX_WARPING_FACTOR,
                                            maxRelXShift = Parameters.MAX_REL_X_SHIFT,
                                            minRelIntervalLength = Parameters.MIN_REL_INTERVAL_LENGTH,
                                            minRelMinMaxDist = Parameters.MIN_REL_MIN_MAX_DIST,
                                            minRelSlopeHeight = Parameters.MIN_REL_SLOPE_HEIGHT,
                                            reference = Parameters.REFERENCE,
                                            outSlope = Parameters.OUT_SLOPE);
        
        # interpolatedXy <- interpolateCurves(consensusProfileData$xWarped, profilesForAlignment, 
        #                                     Parameters.RANKING_DIST_SAMPLE);  # interpolate y-values
        # CurvesMiner.__printAverageDistances(profilesForAlignment, interpolatedXy$y);
        
        # apply the computed x-warping to the real y-coordinates and interpolate 
        # then to a certain number of points
        curve <- getMeanCurve(consensusProfileData$xWarped, 
                              profilesForAlignment, 
                              Defaults.INTERPOLATION_SIZE);
        
        consensusProfilesList <- rlist::list.append(consensusProfilesList, curve$y);
        
        # for tests
        if (plot) Plotter.plotConsensusAndCurves(profilesForAlignment, curve, i, yLimits);
        
      } else {  # else if there is only one curve
        consensusProfilesList <- rlist::list.append(consensusProfilesList, 
                                                    unlist(profiles));
      }
      
      numberOfProfilesList <- rlist::list.append(numberOfProfilesList, numProfiles);
      addedToProfileNumberList <- TRUE;
      years <- rlist::list.append(years, i);
    }
    
    if (!addedToProfileNumberList) numberOfProfilesList <- rlist::list.append(numberOfProfilesList, 0);
  }
  
  return(list(years = years, profiles = consensusProfilesList, numberOfProfiles = numberOfProfilesList));
}


#########################################################
#' Computes the consensus curve.
#'
#' @param characteristicsPerYear {environment("years", "list of profiles")} the profiles per year
#' @param startYear {numerical} the year in which the master-chronology should start
#' @param endYear {numerical} the year in which the master-chronology should end
#'
#' @return {list(years, widths)} the unique years and the widths
#' @export
CurvesMiner.computeCharacteristicsChronology <- function(characteristicsPerYear, startYear, endYear) {
  averagesList <- list();
  yearsList <- list();
  
  # iterate over each year
  for (i in startYear:endYear) {
    characteristics <- characteristicsPerYear[[toString(i)]];
    
    if (!is.null(characteristics)) {  # if for the year width exist
      robustAverage <- Math.tukeysBiweightRobustMean(unlist(characteristics));
      averagesList <- rlist::list.append(averagesList, robustAverage);
      yearsList <- rlist::list.append(yearsList, i);
    }
  }
  
  return(list(years = yearsList, widths = averagesList));
}


#########################################################
#' Plots the mean-distance before and after the alignment.
#'
#' @param beforeAlignment {dataframe} the curves before alignment
#' @param afterAlignment {dataframe} the curves after alignment
CurvesMiner.__printAverageDistances <- function(beforeAlignment, afterAlignment) {
  colnames(afterAlignment) <- 1:ncol(afterAlignment);
  afterAlignment <- as.data.frame.matrix(afterAlignment);
  
  beforeAlignmentMatrix <- CurvesMiner.__getMatrix(beforeAlignment);
  afterAlignmentMatrix <- CurvesMiner.__getMatrix(afterAlignment);
  
  numDistancesBefore <- nrow(beforeAlignmentMatrix);
  numDistancesAfter <- nrow(afterAlignmentMatrix);
  
  diag(beforeAlignmentMatrix) <- 0;
  diag(afterAlignmentMatrix) <- 0;
  
  beforeAlignmentMatrix[lower.tri(beforeAlignmentMatrix)] <- t(beforeAlignmentMatrix)[lower.tri(beforeAlignmentMatrix)];
  afterAlignmentMatrix[lower.tri(afterAlignmentMatrix)] <- t(afterAlignmentMatrix)[lower.tri(afterAlignmentMatrix)];
  
  beforeAlignmentMatrix <- mapply(beforeAlignmentMatrix, FUN = as.numeric);
  afterAlignmentMatrix <- mapply(afterAlignmentMatrix, FUN = as.numeric);
  
  beforeAlignmentMatrix <- matrix(data = beforeAlignmentMatrix, ncol = numDistancesBefore, nrow = numDistancesBefore);
  afterAlignmentMatrix <- matrix(data = afterAlignmentMatrix, ncol = numDistancesAfter, nrow = numDistancesAfter)
  
  distancesBefore <- colSums(beforeAlignmentMatrix);
  distancesAfter <- colSums(afterAlignmentMatrix);
  
  print(sum(distancesBefore) / numDistancesBefore);
  print(sum(distancesAfter) / numDistancesAfter);
  
  print(Strings.SEPERATION);
}


#########################################################
#' Computes a distance matrix for the columns.
#'
#' @param alignment {dataframe} the curves for which distances should be computed
#'
#' @return {matrix} distance matrix
CurvesMiner.__getMatrix <- function(alignment) {
  allDistancesPerProfile <- list();
  
  # iterate over all profiles in that year
  for (j in 1:ncol(alignment)) {
    profile <- alignment[[j]];
    
    distances <- Analyzer.getPairwiseDistances(alignment[-(1:j)], profile, FALSE);
    allDistancesPerProfile <- rlist::list.append(allDistancesPerProfile, distances);
  }
  
  mat <- matrix(list(), nrow = length(allDistancesPerProfile), ncol = length(allDistancesPerProfile));
  
  # set diagonal values
  diag(mat) <- 0;
  
  # iterate over every outer list
  for (i in 1:length(allDistancesPerProfile)) {
    innerList <- allDistancesPerProfile[[i]];
    
    if (length(innerList) > 0) {
      # iterate over all inner list
      for (j in 1:length(innerList)) {
        value <- innerList[[j]];
        
        mat[[i, j+i]] <- value;  # +1 since first position is zero
        mat[[j+i, i]] <- value;
      }
    }
  }
  
  return(mat);
}


#########################################################
#' Subdivides the samples into smaller samples 
#' to avoid a recomputation with smaller samples.
#' Hint: Works only with the given dataset Ostalb due to the name conversion.
#'
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param divisionFactor {numerical} the factor by which the sample length 
#' should be divided to create smaller samples
#' 
#' @return {list} the subdivided sampleData in the same format as before
#' @export
CurvesMiner.subdivideSamples <- function(sampleData, divisionFactor) {
  subdividedSampleData <- list();
  
  newYears <- list(); 
  newProfiles <- list()
  
  # iterate over all samples
  for (i in 1:length(sampleData)) {
    # get the years, profiles and filenames
    years <- sampleData[[i]]$years;
    profiles <- sampleData[[i]]$profiles;
    name <- sampleData[[i]]$name;
    
    # divide year- and sample-lengths by the factor
    numProfiles <- length(profiles);
    partsLength <- floor(numProfiles / divisionFactor);
    
    nameWithoutYears <- stringr::str_sub(name, 1, nchar(name)-Parameters.LENGTH_TO_REMOVE);
    
    for (j in 1:divisionFactor) {
      newYears <- years[(j * partsLength - partsLength + 1):(j * partsLength)];
      newProfiles <- profiles[(j * partsLength - partsLength + 1):(j * partsLength)];
      newName <- paste(nameWithoutYears, newYears[1], 
                       Symbols.HYPHEN, newYears[length(newYears)], sep = Symbols.EMPTY);
      
      sample <- list(years = newYears, profiles = newProfiles, name = newName);
      subdividedSampleData <- rlist::list.append(subdividedSampleData, sample);
    }
  }
  
  return(subdividedSampleData);
}


#########################################################
#' Reconverts loaded samples such that 
#' the CurvesMiner.subtractSampleData-function can be used.
#'
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' Hint: Works only if years at the end of a file name have format ABCD-WXYZ.
#'
#' @return {list(samples, years, names)} the samples with its years and names seperatly
#' @export
CurvesMiner.reconvertSamples <- function(sampleData) {
  # to create
  samples <- list();
  samplesYears <- list();
  curveNames <- list();
  
  # iterate over all samples
  for (i in 1:length(sampleData)) {
    sample <- sampleData[[i]];
    
    samples <- rlist::list.append(samples, sample$profiles);
    samplesYears <- rlist::list.append(samplesYears, sample$years[[1]]);
    curveNames <- rlist::list.append(curveNames, stringr::str_sub(sample$name, 1, stringr::str_length(sample$name)-9));
  }
  
  return(list(samples = samples, years = samplesYears, names = curveNames));
}


#########################################################
#' Extracts corresponding ring widths given the profiles sample data.
#'
#' @param curvesData {list} data of multiple curve-files (ring-widths with its years)
#' @param sampleData {list(samples, years, names)} the profiles with its years and names seperatly
#'
#' @return {list(samples, years, names)} the ring-widths sampleData with its years and names seperatly
#' @export
CurvesMiner.extractCorrespondingWidths <- function(curvesData, sampleData) {
  samples <- list();
  years <- list();
  names <- list();
  
  # iterate overall samples
  for (i in 1:length(sampleData$samples)) {
    sampleLength <- length(sampleData$samples[[i]]);
    sampleYear <- sampleData$years[[i]];
    sampleName <- sampleData$names[[i]];
    
    # search corresponding ring-width in curvesData
    for (j in 1:length(curvesData)) {
      curve <- curvesData[[j]];
      
      if (curve$name == sampleName) {  # if found curve with that name
        
        # iterate over all years in that curve until sample year is found
        for (k in 1:length(curve$years)) {
          year <- curve$years[[k]];
          
          if (sampleYear == year) {  # if year found 
            # extract the widths, years and name
            widths <- curve$widths[k:(k + sampleLength - 1)];
            
            #' year | widths   | name
            #' 1992   w1 w2 w3   "401_MICA-cons"
            #' 1993   w2 w3 w4   "401_MICA-cons"
            samples <- rlist::list.append(samples, widths);
            years <- rlist::list.append(years, sampleYear);
            names <- rlist::list.append(names, sampleName);
          }
        }
      }
    }
  }
  
  return(list(samples = samples, years = years, names = names));
}


#########################################################
#' Subdivides samples into multiple parts which can be tested separately.
#'
#' @param sampleData {list} data of multiple curve-files (profiles or points with its years)
#' @param numDivisions {numeric} the number of divisions
#'
#' @return {list} subdivided data of multiple curve-files (profiles or points with its years)
#' @export
CurvesMiner.subdivideSampleData <- function(sampleData, numDivisions) {
  newSamples <- list();
  
  subsampleLength <- 0;
  
  # iterate over each sample
  for (i in seq_along(sampleData)) {
    sample <- sampleData[[i]];
    
    # extract sample data 
    sampleName <- sample$name;
    sampleProfiles <- sample$profiles;
    sampleYears <- sample$years;
    
    # subdivide sample
    subsampleLength <- length(sampleYears)/numDivisions;
    extractionStart <- 1;
    
    for (j in 1:numDivisions) {
      # extract subsample
      sampleSubProfiles <- sampleProfiles[extractionStart:(extractionStart + subsampleLength - 1)];
      sampleSubYears <- sampleYears[extractionStart:(extractionStart + subsampleLength - 1)];
      
      extractionStart <- extractionStart + subsampleLength;
      
      # create new profile
      newSample <- list();
      
      newSample$name <- sampleName;
      newSample$profiles <- sampleSubProfiles;
      newSample$years <- sampleSubYears;
      
      newSamples <- rlist::list.append(newSamples, newSample);
    }
  }
  
  # create per division a new list
  subdividedSampleData <- list();
  
  # iterate over new created data
  for (i in 1:numDivisions) {
    newSampleData <- list();
    
    for (j in seq_along(sampleData)) {
      newSample <- newSamples[[numDivisions*j + (i-1) - (numDivisions-1)]];
      newSampleData <- rlist::list.append(newSampleData, newSample);
    }
    
    subdividedSampleData <- rlist::list.append(subdividedSampleData, newSampleData);
  }
  
  return(list(subsampleLength = subsampleLength, 
              subdividedSampleData = subdividedSampleData));
}


#########################################################
#' Returns a new data-structure that contains the maximum densities of profiles.
#'
#' @param curvesData {list} the data from multiple curve-files
#'
#' @return {list} the data-structure storing maximum densities of profiles
#' @export
CurvesMiner.getMaximumDensities <- function(curvesData) {
  maximumDensities <- list();
  
  # iterate over all samples
  for (i in seq_along(curvesData)) {
    sample <- curvesData[[i]];
    
    years <- sample$years;
    profiles <- sample$profiles;
    name <- sample$name;
    
    profilesMaxima <- CurvesMiner.__getMaximaOfProfiles(profiles);
    
    densitySample <- list(years = years, maxDensities = profilesMaxima, name = name);
    maximumDensities <- rlist::list.append(maximumDensities, densitySample);
  }
  
  return(maximumDensities);
}


#########################################################
#' Returns the maxima of profiles.
#'
#' @param profiles {list} the profiles for which maxima should be computed
#'
#' @return {list} the maximum-densities per profile
#' @export
CurvesMiner.__getMaximaOfProfiles <- function(profiles) {
  profilesMaxima <- list();
  
  # iterate over all profiles
  for (i in seq_along(profiles)) {
    profile <- unlist(profiles[[i]]);
    maxima <- max(profile);
    profilesMaxima <- rlist::list.append(profilesMaxima, maxima);
  }
  
  return(profilesMaxima);
}