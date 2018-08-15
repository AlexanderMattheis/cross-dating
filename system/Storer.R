#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

STORER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Does the storing of data.

#########################################################
#' Stores the given curve data with its years into one file.
#'
#' @param curveData {list} the curve data of a single tree in all measured directions
#' @param curveFileName {string} the name of the curve
#' @param startYear {numerical} the year in which the data starts
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Stores.storeCurve <- function(curveData, curveFileName, startYear, tableSeparator) {
  allProfiles <- list();
  allYears <- list();
  
  for (i in 1:length(curveData)) {
    profile <- curveData[[i]];
    years <- rep(startYear + (i-1), length(profile));
    
    allProfiles <- rlist::list.append(allProfiles, profile);
    allYears <- rlist::list.append(allYears, years);
  }
  
  data <- data.frame(year = unlist(allYears), consensus = unlist(allProfiles));
  
  write.table(data, file = paste(Paths.OUTPUT, curveFileName, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = tableSeparator);
}


#########################################################
#' Stores samples with consecutive years in an own file to load it later on fastly 
#' or rather avoid recomputation.
#' Hint: Samples are always consecutive. 
#' For broken samples there is another function.
#'
#' @param sampleData {list(samples, years, names)} list of extracted samples with consecutive years
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param widthBased {logical} tells if widths should be used in instead of profiles
#' @export
Storer.storeConsecutiveSamples <- function(sampleData, tableSeparator, widthBased = FALSE) {
  
  # iterate over each sample
  for (i in 1:length(sampleData$samples)) {
    sample <- sampleData$samples[[i]];
    samplelength <- length(sample);
    sampleYear <- as.numeric(sampleData$years[[i]]);
    sampleName <- sampleData$names[[i]];
    
    years <- sampleYear:(sampleYear + samplelength - 1);
    
    curveFileName <- paste(sampleName, years[[1]], 
                           Symbols.HYPHEN, years[[length(years)]], 
                           Extensions.CSV, sep = Symbols.EMPTY);
    
    if (widthBased) {  # extendions for ring widths
      curveData <- list(years = years, widths = sample);
      Storer.storeWidths(curveData, curveFileName, tableSeparator);
    } else {
      curveData <- list(years = years, profiles = sample);
      Storer.storeCurve2(curveData, curveFileName, tableSeparator);
    }
  }
}


#########################################################
#' Stores samples with consecutive years in an own file to load it later on fastly 
#' or rather avoid recomputation.
#'
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeConsecutiveSamples2 <- function(sampleData, tableSeparator) {
  
  # iterate over each sample
  for (i in 1:length(sampleData)) {
    sample <- sampleData[[i]]$profiles;
    years <- as.list(as.numeric(unlist(sampleData[[i]]$years)));
    sampleName <- paste(sampleData[[i]]$name, Extensions.CSV, sep = Symbols.EMPTY);
    
    curveData <- list(years = years, profiles = sample);
    Storer.storeCurve2(curveData, sampleName, tableSeparator);
  }
}


#########################################################
#' Stores each sample in an own file to load it later on fastly 
#' or rather avoid recomputation.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeSamples2 <- function(curvesData, tableSeparator) {
  # iterate over each sample
  for (i in 1:length(curvesData)) {
    profiles <- curvesData[[i]]$profiles;
    years <- curvesData[[i]]$years;
    curveName <- curvesData[[i]]$name;
    
    curveData <- list(years = years, profiles = profiles);
    curveFileName <- paste(curveName, Extensions.CSV, sep = Symbols.EMPTY);
    Storer.storeCurve2(curveData, curveFileName, tableSeparator);
  }
}


#########################################################
#' Stores the given curve data with its years into one file.
#'
#' @param curveData {list(years, profiles)} curve data with years
#' @param curveFileName {string} the name you want give the curve
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param colNames {vector} vector of length two with the names for both columns
#' @export
Storer.storeCurve2 <- function(curveData, curveFileName, tableSeparator, colNames = vector()) {
  allProfiles <- list();
  allYears <- list();
  
  for (i in 1:length(curveData$years)) {
    profile <- curveData$profiles[[i]];
    years <- rep(curveData$years[[i]], length(profile));
    
    allProfiles <- rlist::list.append(allProfiles, profile);
    allYears <- rlist::list.append(allYears, years);
  }
  
  data <- data.frame(year = unlist(allYears), GD = unlist(allProfiles));  # GD: gravimetric density
  
  if (length(colNames) == 2) {  # extension fo interface
    colnames(data) <- colNames;
  }
  
  write.table(data, file = paste(Paths.OUTPUT, curveFileName, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = tableSeparator);
}


#########################################################
#' Stores the widths.
#'
#' @param curveData {list(years, widths)} curve data with years
#' @param curveFileName {string} the name you want give the curve
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeWidths <- function(curveData, curveFileName, tableSeparator) {
  data <- data.frame(year = unlist(curveData$years), width = unlist(curveData$widths));
  
  write.table(data, file = paste(Paths.OUTPUT, curveFileName, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = tableSeparator);
}


#########################################################
#' Stores all computed scores for each sample.
#' Deprecated but needed for compatibility since R conversion issues
#' with dataframes.
#'
#' @param scoresPerSample {list} the scores per sample
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param chronologyYears {list} the years of the master-chronology
#' @param method {string} the method which was used [a, b, c, d, e, f, g, h]
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeScores <- function(scoresPerSample, sampleData, chronologyYears, method, tableSeparator) {
  # iterate over all samples
  for (i in 1:length(scoresPerSample)) {
    # retrieve scores and sample-name
    scores <- scoresPerSample[[i]];
    sampleName <- sampleData[[i]]$name;
    numberOfScores <- length(scores);
    years <- chronologyYears[1:numberOfScores];
    
    # create file with that sample-name
    data <- data.frame(year = unlist(years), score = unlist(scores)); 
    
    write.table(data, file = paste(Paths.OUTPUT, sampleName, Strings.SCORE, Symbols.UNDERSCORE,
                                   method, Extensions.CSV, sep = Symbols.EMPTY), 
                row.names = FALSE, sep = tableSeparator);
  }
}


#########################################################
#' Stores individual scores for each sample.
#'
#' @param scoresPerSample {list} the scores per sample
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param chronologyYears {list} the years of the master-chronology
#' @param method {string} the method which was used [a, b, c, d, e, f, g, h]
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeIndividualScores <- function(scoresPerSample, sampleData, chronologyYears, method, tableSeparator) {
  # iterate over all samples
  for (i in 1:length(scoresPerSample)) {
    # retrieve scores and sample-name
    scores <- scoresPerSample[[i]];
    sampleName <- sampleData[[i]]$name;
    numberOfScores <- length(scores);
    years <- chronologyYears[1:numberOfScores];
    
    # create file with that sample-name
    data <- do.call(rbind, scores);
    data <- cbind(years, data);
    
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
    
    write.table(data, file = paste(Paths.OUTPUT, sampleName, Strings.SCORE, Symbols.UNDERSCORE,
                                   method, Extensions.CSV, sep = Symbols.EMPTY), 
                row.names = FALSE, sep = tableSeparator);
  }
}


#########################################################
#' Stores all values from a list.
#'
#' @param list {list} the values which you want store
#' @param fileName {string} the name you want give the file
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @param quotes {logical} tells if strings should be surrounded with quotes
#' @export
Storer.storeList <- function(list, fileName, tableSeparator, quotes = TRUE) {
  data <- data.frame(list = unlist(list)); 
  write.table(data, file = paste(Paths.OUTPUT, fileName, Extensions.CSV, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = tableSeparator, quote = quotes)
}


#########################################################
#' Stores all values from a list of lists as a matrix.
#' Only the row names are stored, since they are identical to the column names.
#'
#' @param listOfList {list} the list which stores a list
#' @param names {list} the strings for the columns or rows of the matrix
#' @param fileName {string} the name of the file to store
#' @export
Storer.storeMatrix <- function(listOfList, names, fileName) {
  mat <- matrix(list(), nrow = length(names), ncol = length(names));
  
  # set diagonal values
  diag(mat) <- 0;
  
  # iterate over every outer list
  for (i in 1:length(listOfList)) {
    innerList <- listOfList[[i]];
    
    if (length(innerList) > 0) {
      # iterate over all inner list
      for (j in 1:length(innerList)) {
        value <- innerList[[j]];
        
        mat[[i, j+i]] <- value;  # +1 since first position is zero
        mat[[j+i, i]] <- value;
      }
    }
  }
  
  colnames(mat) <- unlist(names);
  write.table(mat, file = paste(Paths.OUTPUT, fileName, sep = Symbols.EMPTY), 
              sep = Symbols.REAL_DATA_SEPARATOR, row.names = FALSE);
}


#########################################################
#' Stores the plot on the disk.
#'
#' @param plot {ggplot} the plot you want save on hard drive
#' @param filename {string} the name you want give the file
#' @export
Storer.savePlot <- function(plot, filename) {
  ggsave(filename = filename, path = Paths.OUTPUT,  plot = plot);
}


#########################################################
#' Stores the encoded plot on the disk.
#'
#' @param plot {ggplot} the plot you want save on hard drive
#' @param patternLength {numerical} the current pattern length for which scores are plotted
#' @param parametersToEncode {list} the parameters which shoulde be encoded into plot
Storer.saveEncodedPlot <- function(plot, patternLength, parametersToEncode) {
  distanceFunction <- parametersToEncode$distanceFunction;
  loadedPattern <- parametersToEncode$loadedPattern;
  subpatternStart <- parametersToEncode$subpatternStart;
  subpatternEnd <- parametersToEncode$subpatternEnd;
  treeForConsensus <- parametersToEncode$treeForConsensus;
  
  filename <- paste(loadedPattern, distanceFunction, subpatternStart, subpatternEnd, treeForConsensus, patternLength);
  
  ggsave(filename = paste(filename, Extensions.PDF, sep=Symbols.EMPTY), 
         path = Paths.OUTPUT, plot = plot, device = "pdf");
}


#########################################################
#' Stores the Venn diagram on the hard disk together with other data.
#'
#' @param maximumIntesection {list} the intersection of all four methods 
#' @param pass {numerical} the pass to encode into file names 
#' @export
Storer.saveVennDiagram <- function(maximumIntesection, pass) {
  filename <- paste(Files.VENN_DIAGRAM, pass, Extensions.TIFF, sep = Symbols.EMPTY);
  
  maximumIntesection <- gsub(Strings.TO_REPLACE_1, Symbols.EMPTY, as.list(maximumIntesection));
  maximumIntesection <- gsub(Strings.TO_REPLACE_2, Strings.REPLACEMENT, as.list(maximumIntesection));
  
  Storer.storeList(maximumIntesection, 
                   paste(Files.INTERSECTION, pass, sep = Symbols.EMPTY), 
                   Symbols.REAL_DATA_SEPARATOR, quotes = FALSE);
  tiff(filename = filename, compression = "lzw");
  grid.draw(venn.plot);
  dev.off();
}


#########################################################
#' Inserts characteristics like ring-widths and corresponding profiles into a single file.
#' 
#' @param pathProfiles {string} the path to the profiles
#' @param pathCharacteristics {string} the path to the widths/maximum densities
#' @param profileName {string} the name of the file containing profiles (without extension)
#' @param widthName {string} the name of the file containing widths (without extension)
#' @param outputFileName {string} the name for the output (without extension)
#' @param extension {string} the extension of the files
#' @param isSample {logical} tells if it is a sample
#'
#' @examples
#' #consensus of profiles [GIVEN]
#' "year","GD"
#' 1992, 0.23
#' ...
#' 1992, 0.54
#' 1993, 0.17
#' ...
#' 1993, 0.67
#' ...
#' 
#' #consensus of characteristic [GIVEN]
#' "year","width"
#' 1992, 0.35
#' 1993, 0.52
#' ...
#' 
#' #consensus [DEFAULT]
#' year, density, characteristic (optional)
#' 1992, 0.23,    0.35
#' ...
#' 1992, 0.54,    0.35
#' 1993, 0.17,    0.52
#' ...
#' 1993, 0.67,    0.52
#' ...
#' 
#' #samples [DEFAULT]
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
Storer.storeInDefaultFormat <- function(pathProfiles, pathCharacteristics, 
                                        profileName, widthName, outputFileName, 
                                        extension = Extensions.CSV, isSample = FALSE) {
  
  # {list(years, profiles, name)} 
  profilesData <- Loader.readInProfiles(pathProfiles, profileName, Symbols.REAL_DATA_SEPARATOR, extension);
  
  # {list(years, widths, name)} 
  characteristicData <- Loader.readInWidths(pathCharacteristics, widthName, Symbols.REAL_DATA_SEPARATOR, extension);
  
  # retrieve data
  name <- characteristicData$name;
  
  years <- profilesData$years;
  profiles <- profilesData$profiles;
  characteristicValues <- characteristicData$widths;
  
  # iterate over all years and create [DEFAULT] consensus
  newYears <- list();
  newProfiles <- list();
  newWidths <- list();
  
  for (i in 1:length(years)) {
    year <- as.numeric(years[[i]]);
    profile <- profiles[[i]];
    characteristic <- as.numeric(characteristicValues[[i]]);
    
    profileLength <- length(profile);
    
    newYears <- rlist::list.append(newYears, rep(year, profileLength));
    newProfiles <- rlist::list.append(newProfiles, profile);
    newWidths <- rlist::list.append(newWidths, rep(characteristic, profileLength));                                  
  }

  data <- data.frame();
  
  if (isSample) {  # convert to parts (sample store parts since the years are unknown)
    newYearsVec <- unlist(newYears);
    parts <- newYearsVec - min(newYearsVec) + 1;
    data <- data.frame(part = parts, density = unlist(newProfiles), characteristic = unlist(newWidths)); 
  } else {
    data <- data.frame(year = unlist(newYears), density = unlist(newProfiles), characteristic = unlist(newWidths)); 
  }
  
  write.table(data, file = paste(Paths.OUTPUT, outputFileName, extension, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Inserts ring widths and corresponding profiles into a single file.
#' Hint: Files in both folders must have same names!
#' 
#' @param pathProfiles {string} the path to the profiles
#' @param pathCharacteristics {string} the path to the widths/maximum densities
#' @param isSample {logical} tells if it is a sample
#' @export
Storer.storeAllInDefaultFormat <- function(pathProfiles, pathCharacteristics, isSample = TRUE) {
  fileNames <- list.files(path = pathProfiles);

  for (i in 1:length(fileNames)) {
    fileName <- fileNames[[i]];
    
    Storer.storeInDefaultFormat(pathProfiles, pathCharacteristics, 
                                fileName, fileName, fileName, Symbols.EMPTY, isSample);
  }
}


#########################################################
#' Stores the given matrix on the harddrive.
#'
#' @param table {matrix} the matrix you want store
#' @param fileName {string} the name of the file to store
#' @param extension {string} the extension which should be used
#' @export
Storer.storeTable <- function(table, fileName, extension) {
  write.table(table, file = paste(Paths.OUTPUT, fileName, extension, sep = Symbols.EMPTY), 
              sep = Symbols.REAL_DATA_SEPARATOR, row.names = FALSE, quote = FALSE);
}


#########################################################
#' Removes from the consensus the given years and stores in default format.
#'
#' @param consensusPath {string} the path to the consensus
#' @param consensusName {string} the name of the consensus file
#' @param yearsToRemove {vector} the years to remove
#' @param outputFileName {string} the name for the output (without extension)
#' @export
Storer.removeYears <- function(consensusPath, consensusName, yearsToRemove, outputFileName) {
  # list(years = curveYears, profiles = profiles, name = curveName)
  consensus <- Loader.readInDefaultConsensus(consensusPath, consensusName, 
                                             Symbols.REAL_DATA_SEPARATOR, 
                                             Extensions.CSV);
  
  # remove if there is something to remove
  if (length(yearsToRemove) > 0) {
    toRemove <- -which(consensus$years %in% yearsToRemove);
    
    consensus$years <- consensus$years[toRemove];  # do not change the order!
    consensus$profiles <- consensus$profiles[toRemove];
    consensus$characteristics <- consensus$characteristics[toRemove];
  } 
  
  # encode
  allYears <- list();
  allProfiles <- list();
  allWidths <- list();
  
  for (i in 1:length(consensus$years)) {
    profile <- consensus$profiles[[i]];
    years <- rep(consensus$years[[i]], length(profile));
    characteristics <-  rep(consensus$characteristics[[i]], length(profile));
    
    allYears <- rlist::list.append(allYears, years);
    allProfiles <- rlist::list.append(allProfiles, profile);
    allWidths <- rlist::list.append(allWidths, characteristics);
  }
  
  data <- data.frame(year = unlist(allYears), density = unlist(allProfiles), characteristic = unlist(allWidths));
  
  write.table(data, file = paste(Paths.OUTPUT, outputFileName, Extensions.CSV, sep = Symbols.EMPTY), 
              row.names = FALSE, sep = Symbols.REAL_DATA_SEPARATOR);
}


#########################################################
#' Removes from the consensi the given samples and stores in default format.
#'
#' @param consensiPath {string} the path to the per tree consensi
#' @param samplesPath {string} the path to the samples which have to be deleted from the consensi
#' @export
Storer.removeSampleYearsFromAll <- function(consensiPath, samplesPath) {
  consensiNames <- list.files(path = consensiPath);
  
  # iterate oer all consensi
  for (i in 1:length(consensiNames)) {
    consensiName <- consensiNames[[i]];
    consensiName <- stringr::str_sub(consensiName, start = 1, end = nchar(consensiName)-nchar(Extensions.CSV));
    patternName <- paste(Expression.START, consensiName, sep = Symbols.EMPTY);  # to avoid selecting with pattern 401_... the pattern 1401_..
    
    # search in samples path for that consensi name
    samplesNames <- list.files(path = samplesPath, pattern = patternName);
    
    yearsToRemove <- list();
    
    if (length(samplesNames) > 0) {
      # iterate over all that sample names
      for (j in 1:length(samplesNames)) {
        sampleName <- samplesNames[[j]];
        years <- stringr::str_extract(sampleName, Expression.YEARS_IN_BETWEEN);
        yearsSplitted <- unlist(stringr::str_split(years, Symbols.HYPHEN));
        
        startYear <- yearsSplitted[[1]];
        endYear <- yearsSplitted[[2]];
        
        yearsToRemove <- rlist::list.append(yearsToRemove, startYear:endYear);
      }
    }

    Storer.removeYears(consensusPath = consensiPath, 
                       consensusName = consensiName, 
                       yearsToRemove = unlist(yearsToRemove),
                       outputFileName = consensiName);
  }
}


#########################################################
#' Stores the profiles in the given environment.
#'
#' @param environment {environment} the environment to store
#' @param keys {vector} the keys to store
#' @export
Storer.storeProfilesEnvironment <- function(environment, keys) {
  # iterate over each key that should be stored
  for (key in keys) {
    profiles <- environment[[toString(key)]];
    numericProfiles <- list();
    
    for (i in seq_len(length(profiles))) {
      numericProfile <- unlist(profiles[[i]]);
      numericProfiles <- rlist::list.append(numericProfiles, numericProfile);
    }
    
    data <- do.call(cbind.data.frame, numericProfiles);
    
    write.table(data, file = paste(Paths.OUTPUT, key, Extensions.CSV, sep = Symbols.EMPTY), 
                row.names = FALSE, col.names = FALSE, sep = Symbols.REAL_DATA_SEPARATOR);
  }
}


#########################################################
#' Stores the counts per year of each sample.
#'
#' @param countsPerSample {list} the counts per sample
#' @param sampleData {list} data of multiple curve-files (profiles with its years)
#' @param chronologyYears {list} the years of the master-chronology
#' @param method {string} the method which was used [a, b, c, d]
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeCounts <- function(countsPerSample, sampleData, chronologyYears, method, tableSeparator) {
  # iterate over all samples
  for (i in 1:length(countsPerSample)) {
    # retrieve scores and sample-name
    counts <- countsPerSample[[i]];
    sampleName <- sampleData[[i]]$name;
    numberOfCounts <- length(counts);
    years <- chronologyYears[1:numberOfCounts];
    
    # create file with that sample-name
    data <- do.call(rbind, counts);
    data <- cbind(years, data);
    
    colnames(data) <- c(Strings.YEAR_LOWER_CASE, Strings.COUNT_LOWER_CASE);
    
    write.table(data, file = paste(Paths.OUTPUT, sampleName, Strings.COUNT_2, Symbols.UNDERSCORE,
                                   method, Extensions.CSV, sep = Symbols.EMPTY), 
                row.names = FALSE, sep = tableSeparator);
  }
}


#########################################################
#' Stores the maximum densities of profiles.
#'
#' @param maximumDensities {list} the maximum densities of profiles
#' @param tableSeparator {string} the separator-symbol in the curve csv data
#' @export
Storer.storeMaxDensities <- function(maximumDensities, tableSeparator) {
  # iterate over all samples
  for (i in 1:length(maximumDensities)) {
    sample <- maximumDensities[[i]];
    
    # read out data
    years <- sample$years;
    maxDensities <- sample$maxDensities;
    sampleName <- sample$name;
    
    # create file with that sample-name
    data <- do.call(rbind, maxDensities);
    data <- cbind(years, data);
    
    colnames(data) <- c(Strings.YEAR_LOWER_CASE, Defaults.DENSITY_NAME);
    
    write.table(data, file = paste(Paths.OUTPUT, sampleName, Extensions.CSV, sep = Symbols.EMPTY),
                row.names = FALSE, sep = tableSeparator, quote = FALSE)
  }
}