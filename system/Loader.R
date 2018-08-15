#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

LOADER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Does the loading of curve data and scores.

#########################################################
#' Given curves by a path, the individual CSV-curves with equidistant x-coordinates are read in.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param readWidths {logical} tells if widths should be read in instead of profiles
#' @param extension {string} the extension which should be used for the files
#' @param columnName {string} the name of the column with the profile-data
#' @param columnNameWidth {string} the name of the column with the ring-width-data or maximum density
#' 
#' @return {list} data of multiple curve-files (profiles or points with its years)
#' @export
Loader.readInCurves <- function(datasetPath, tableSeparator, readWidths = FALSE, extension = Extensions.CSV,
                                columnName = Defaults.GENERIC_DENSITY_COLUMN_NAME, 
                                columnNameWidth = Defaults.GENERIC_WIDTH_COLUMN_NAME) {
  paths <- list.files(path = datasetPath);
  # View(paths);
  curvesData <- list();
  data <- list();
  
  for (i in 1:length(paths)) {
    fullName <- paths[[i]];  # with extension
    
    # check-ups
    if (!Loader.correctExtension(fullName, extension)) stop(Exceptions.WRONG_EXTENSION);
    
    curveName <- substr(fullName, 1, nchar(fullName)-nchar(extension));  # without extension
    
    if (readWidths) {  # stored in separate file
      data <- Loader.readInWidths(datasetPath, curveName, tableSeparator, extension,
                                  columnName = columnNameWidth);
    } else {
      data <- Loader.readInProfiles(datasetPath, curveName, tableSeparator, extension, 
                                    columnName = columnName);
    }
    
    curvesData <- rlist::list.append(curvesData, data);
  }
  
  return(curvesData);
}


#########################################################
#' Checks the extension.
#'
#' @param fullName {string} the path to the file
#' @param desiredExtension {string} the desired extension
#'
#' @return {logical} true if the extension of the file equals the desired extension
#' @export
Loader.correctExtension <- function(fullName, desiredExtension) {
  extension <- substr(fullName, nchar(fullName)-nchar(desiredExtension) + 1, nchar(fullName));
  return(identical(extension, desiredExtension));
}


#########################################################
#' Given a curve by a path, the individual 
#' year widths are extracted and returned.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param curveName {string} the curve from which profiles should be extracted
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param extension {string} the extension which should be used for ring width curves
#' @param columnName {string} the name of the characteristic column
#'
#' @return {list(years, widths, name)} the unique years and the profiles
#' @export
Loader.readInWidths <- function(datasetPath, curveName, tableSeparator, extension = Extensions.META, 
                                columnName = Defaults.GENERIC_WIDTH_COLUMN_NAME) {
  # read in data
  path <- paste(datasetPath, curveName, extension, sep = Symbols.EMPTY);
  curveData <- na.omit(read.csv(path, header = TRUE, sep = tableSeparator));  # omit rows with NAs
  
  # extract
  if (length(curveData$meanNumPoints) == 0) {
    widths <- as.list(curveData[[columnName]]);
  } else {
    widths <- as.list(curveData$meanNumPoints);
  }
  years <- as.list(curveData$year);

  return(list(years = years, widths = widths, name = curveName));
}


#########################################################
#' Given a curve by a path, the individual 
#' year profiles are extracted and returned.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param curveName {string} the curve from which profiles should be extracted
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param extension {string} the extension which should be used for profile curves
#' @param columnName {string} the name of the characteristic column
#' 
#' @return {list(years, profiles, name)} the unique years and the profiles
#' @export
Loader.readInProfiles <- function(datasetPath, curveName, tableSeparator, 
                                  extension = Extensions.CSV, columnName = Defaults.GENERIC_DENSITY_COLUMN_NAME) {
  # read in data
  path <- paste(datasetPath, curveName, extension, sep = Symbols.EMPTY);
  curveData <- na.omit(read.csv(path, header = TRUE, sep = tableSeparator));  # omit rows with NAs
  curveLength <- length(curveData$year);
  
  # check equidistance
  # relXinRing <- curveData$relXinRing;
  # equidistant <- TRUE;
  
  curveYears <- levels(factor(curveData$year));
  numberOfYears <- length(curveYears);  # hint: removed years possible 

  j <- 1;  # vertical position in table
  
  # check steadiness in years
  lastYear <- 0;
  if (numberOfYears > 0) lastYear <- as.numeric(curveData$year[[1]]);
  steady <- TRUE;
  
  # data to store
  profiles <- list();
  
  # iterate overall all years

  for (i in 1:numberOfYears) {
    
    profile <- list();
    currentYear <- as.numeric(curveYears[[i]]);
    
    # while currentYear == curveData$year[[j]] add values to profile
    while(j <= curveLength && currentYear == curveData$year[[j]]) {
      if (lastYear > as.numeric(curveData$year[[j]])) steady <- FALSE;
      
      if (is.null(curveData[[curveName]])) {  # in case of real data
        profile <- rlist::list.append(profile, 
                                      curveData[[columnName]][[j]]);
      } else {  # in case of artificial data
        profile <- rlist::list.append(profile, curveData[[curveName]][[j]]);
      }
      
      lastYear <- as.numeric(curveData$year[[j]]);
      j = j + 1; 
    }
    
    positions <- curveLength - 1;
    # if (j < curveLength && round(relXinRing[j:(j + curveLength - 1)], 4) != round((0:positions)/positions, 4)) {
    #   equidistant <- FALSE;
    # }
    
    profiles <- rlist::list.append(profiles, profile);
  }

  # print(Titles.STEADY);
  # print(steady);
  # print(Strings.SEPERATION);
  # print(Titles.EQUIDISTANT);
  # print(equidistant);
  # print(Strings.SEPERATION);
  # print(Titles.PROFILES_AND_YEARS);
  # print(length(relXinRing));
  # print(numberOfYears);
  # print(Strings.SEPERATION);
  # print(Titles.CONSECUTIVE);
  # print(curveName);
  # print(length(curveYears) == (as.numeric(curveYears[[length(curveYears)]]) - as.numeric(curveYears[[1]]) + 1 ));
  # print(Strings.SEPERATION);
  
  return(list(years = curveYears, profiles = profiles, name = curveName));
}


#########################################################
#' Reads in score files with a single score-column per file. 
#'
#' @param scoresSetPath {string} the path to a set of scores
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param fullNames {logical} tells if the full filenames with the full path should be returned or not
#'
#' @return {list(scoresData, fileNames, scoreTables)} per sample scores and filenames 
#' (also returning the per sample scores table)
#' @export
Loader.readInScores <- function(scoresSetPath, tableSeparator, fullNames = FALSE) {
  paths <- list.files(path = scoresSetPath, full.names = TRUE);  # full.names = TRUE: prepend path
  # View(paths);
  
  scoresData <- list();
  fileNames  <- list();
  scoresTables <- list();
  
  for (i in 1:length(paths)) {
    filePath <- paths[[i]];  # with extension and prepended path
    #print(filePath);
    
    if (fullNames) {
      fullName <- filePath;
    } else {  # get name an year
      tokens <- stringr::str_split(filePath, Paths.PATH)[[1]];
      fullName <- tokens[[length(tokens)]];
    }
    
    fileNames <- rlist::list.append(fileNames, fullName);
    
    # save stored data
    entries <- read.csv(filePath, header = TRUE, sep = tableSeparator);
    scoresTables <- rlist::list.append(scoresTables, entries);
    scoresData <- rlist::list.append(scoresData, entries$score);
  }
  
  return(list(scoresData = scoresData, fileNames = fileNames, scoresTables = scoresTables));
}


#########################################################
#' Reads individual score-columns.
#'
#' @param scoresSetPath {string} the path to a set of scores
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#'
#' @return {list(scoresTables, fileNames)} per sample scores and file names
#' @export
Loader.readInIndividualScores <- function(scoresSetPath, tableSeparator) {
  paths <- list.files(path = scoresSetPath, full.names = TRUE);  # full.names = TRUE: prepend path
  #View(paths);
  
  scoresTables <- list();
  fileNames  <- list();
  
  for (i in 1:length(paths)) {
    filePath <- paths[[i]];  # with extension and prepended path
    #print(filePath);
    
    # get name an year
    tokens <- stringr::str_split(filePath, Paths.PATH)[[1]];
    fullName <- tokens[[length(tokens)]];
    
    fileNames <- rlist::list.append(fileNames, fullName);
    
    # save stored data
    entries <- read.csv(filePath, header = TRUE, sep = tableSeparator);
    scoresTables <- rlist::list.append(scoresTables, entries);
  }
  
  return(list(scoresTables = scoresTables, fileNames = fileNames));
}


#########################################################
#' Reads in list from a given file. 
#'
#' @param filePath {string} the path to the file
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#'
#' @return {list} per year data
Loader.readInList <- function(filePath, tableSeparator) {
  entries <- read.csv(filePath, header = TRUE, sep = tableSeparator);
  return(entries$list);
}


#########################################################
#' Reads in matrix from a given file.
#'
#' @param filePath {string} the path to the file
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#'
#' @return {dataframe} the matrix as an dataframe
Loader.readInMatrix <- function(filePath, tableSeparator) {
  entries <- read.csv(filePath, header = TRUE, sep = tableSeparator);
  return(entries);
}


#########################################################
#' Reads in sample files.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param extension {string} the extension which should be used for the files
#'
#' @return {list} data of multiple curve-files (samples with all its data)
#' @example files
#' #samples
#' part,    density,  characteristic
#' 1,       0.17,     0.51
#' 1,       0.23,     0.51
#' ...      ... ,     ...
#' 1,       0.57,     0.51
#' 2,       0.18,     0.65
#' 2,       0.19,     0.65
#' ...      ... ,     ...
#' 2,       0.59,     0.65
#' ...
#' 
#' @export
Loader.readInSamples <- function(datasetPath, tableSeparator, extension) {
  paths <- list.files(path = datasetPath);
  # View(paths);
  curvesData <- list();
  data <- list();
  
  for (i in 1:length(paths)) {
    fullName <- paths[[i]];  # with extension
    
    # check-up
    if (!Loader.correctExtension(fullName, extension)) stop(Exceptions.WRONG_EXTENSION);
    
    # load
    curveName <- substr(fullName, 1, nchar(fullName)-nchar(extension));  # without extension
    data <- Loader.readInSample(datasetPath, curveName, tableSeparator, extension);

    # append
    curvesData <- rlist::list.append(curvesData, data);
  }
  
  return(curvesData);
}


#########################################################
#' Reads in sample file.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param curveName {string} the curve from which profiles should be extracted
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param extension {string} the extension which should be used for the files
#'
#' @return {list(parts, profiles, widths, name)} the unique years and the profiles and years
#' @example file
#' #samples
#' part,    density,  characteristic
#' 1,       0.17,     0.51
#' 1,       0.23,     0.51
#' ...      ... ,     ...
#' 1,       0.57,     0.51
#' 2,       0.18,     0.65
#' 2,       0.19,     0.65
#' ...      ... ,     ...
#' 2,       0.59,     0.65
#' ...
#' 
#' @export
Loader.readInSample <- function(datasetPath, curveName, tableSeparator, extension, partColumnName = Defaults.PART_NAME) {
  # read in data
  path <- paste(datasetPath, curveName, extension, sep = Symbols.EMPTY);
  curveData <- na.omit(read.csv(path, header = TRUE, sep = tableSeparator));  # omit rows with NAs
  curveLength <- length(curveData[[partColumnName]]);
  
  curveParts <- levels(factor(curveData[[partColumnName]]));
  curveCharacteristics <- curveData$characteristic;
  numberOfParts <- length(curveParts);
  
  j <- 1;  # vertical position in table

  # data to store
  parts <- list();
  profiles <- list();
  characteristics <- list();
  
  # iterate overall all years
  for (i in 1:numberOfParts) {
    profile <- list();
    currentPart <- as.numeric(curveParts[[i]]);
    
    if(!is.null(curveCharacteristics)) {  # since widths are optional
      currentCharacteristics <- as.numeric(curveCharacteristics[[j]]);
      characteristics <- rlist::list.append(characteristics, currentCharacteristics);
    }
    
    # read in profile: while currentPart == curveData$part[[j]] add value to profile
    while(j <= curveLength && currentPart == curveData[[partColumnName]][[j]]) {
      profile <- rlist::list.append(profile, 
                                    curveData[[Defaults.DENSITY_NAME]][[j]]);
      j = j + 1; 
    }
    
    parts <- rlist::list.append(parts, currentPart);
    profiles <- rlist::list.append(profiles, profile);
  }
  
  return(list(parts = parts, profiles = profiles, characteristics = characteristics, name = curveName));
}


#########################################################
#' Given a curve by a path, the individual 
#' year profiles and widths are extracted and returned.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param curveName {string} the curve from which profiles should be extracted
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param extension {string} the extension which should be used for profile curves
#' 
#' @return {list(years, profiles, widths name)} the unique years and the profiles together with the corresponding widths
#' @example 
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
#' @export
Loader.readInDefaultConsensus <- function(datasetPath, curveName, 
                                          tableSeparator, extension) {
  data <- Loader.readInSample(datasetPath, curveName, 
                              tableSeparator, extension, 
                              partColumnName = Defaults.YEAR_NAME);
  
  years <- data$parts;
  profiles <- data$profiles;
  characteristics <- data$characteristics;
  
  return(list(years = years, profiles = profiles, characteristics = characteristics, name = curveName));
}


#########################################################
#' Given default consensi by a path, 
#' the individual CSV-curves are read in.
#'
#' @param datasetPath {string} the path to a set of curves
#'
#' @return {list} data of multiple curve-files (profiles and ring-widths with its years)
#' @example 
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
#' @export
Loader.readInDefaultConsensi <- function(datasetPath) {
  paths <- list.files(path = datasetPath);
  # View(paths);
  curvesData <- list();
  data <- list();
  
  for (i in 1:length(paths)) {
    fullName <- paths[[i]];  # with extension
    
    data <- Loader.readInDefaultConsensus(datasetPath = datasetPath, 
                                          curveName = fullName, 
                                          tableSeparator = Symbols.REAL_DATA_SEPARATOR, 
                                          extension = Symbols.EMPTY);
    
    curvesData <- rlist::list.append(curvesData, data);
  }
  
  return(curvesData);
}


#########################################################
#' Loads files containing multiple profiles as rows.
#'
#' @param datasetPath {string} the path to a set of curves
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#'
#' @return {list(perYearProfiles, startYear, endYear)} the environment containing profiles per year
#' @export
Loader.loadBuckets <- function(datasetPath, tableSeparator) {
  paths <- list.files(path = datasetPath, full.names = TRUE);
  
  perYearProfiles <- new.env();
  startYear <- -1;
  endYear <- -1;
  
  for (i in 1:length(paths)) {
    path <- paths[[i]];
    splits <- unlist(stringr::str_split(path, Symbols.SLASH));
    lastSplit <- splits[[length(splits)]];
    name <- substr(lastSplit, 1, nchar(lastSplit)-nchar(Extensions.CSV));
    
    if (i == 1) {
      startYear <- as.numeric(name);
    } else if (i == length(paths)) {
      endYear <- as.numeric(name);
    }
    
    profiles <- as.list(read.csv(path, header = FALSE, sep = tableSeparator));
    profiles <- lapply(profiles, as.list);
    perYearProfiles[[name]] <- profiles;
  }
  
  return(list(profilesPerYear = perYearProfiles, minStartYear = startYear, maxStartYear = endYear));
}


#########################################################
#' Reads in counts files. 
#'
#' @param countsPath {string} the path to a set of per year counts
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param fullNames {logical} tells if the full filenames with the full path should be returned or not
#'
#' @return {list} the counts (per sample)
#' @export
Loader.readInCounts <- function(countsPath, tableSeparator, fullNames = FALSE) {
  paths <- list.files(path = countsPath, full.names = TRUE);  # full.names = TRUE: prepend path
  
  countsData <- list();
  
  for (i in 1:length(paths)) {
    filePath <- paths[[i]];  # with extension and prepended path
    
    # store 
    entries <- read.csv(filePath, header = TRUE, sep = tableSeparator);
    countsData <- rlist::list.append(countsData, entries$count);
  }
  
  return(countsData);
}