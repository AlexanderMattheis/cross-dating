#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

DATA_ANALYSIS_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Executes the data analysis which was used to create the Data Acquisition chapter.

#########################################################
#' Starts the analysis.
DataAnalysis.start <- function() {
  # check if for each point exactly one profile exists 
  # by checking if the histograms are equal (points per year == profiles per year) - just for safety
  # DataAnalysis.__checkPointsPerYear(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR);
  
  # create histogram of dataset
  # DataAnalysis.__createProfilesPerYear(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR);
  
  # store distances within data
  # DataAnalysis.__computePairwiseDistancesPerYear(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR);
  
  # plot distances
  # DataAnalysis.__createAllVsAllDistancesPlot(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR);
  # DataAnalysis.__createAllVsAllDistancesPlot(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR, TRUE,
  #                                          Titles.PROFILES_MEAN_DISTANCES_PER_YEAR);
  
  # compute consenus-curve plot
  # DataAnalysis.__createConsensus(Defaults.CURVE_OSTALB_START_YEAR, Defaults.CURVE_OSTALB_END_YEAR);
  
  # compare neighboured consensi profiles
  # DataAnalysis.__compareCorrelationCoefficients(Paths.INPUT, Files.OSTALB_CONSENSUS);
}


#########################################################
#' Checks if the number of points in a year is really equal 
#' to the number of profiles in the year (just for safety).
#'
#' @param startYear {number} the year at which the check up should start
#' @param endYear {number} the year at which the check up should end
#' @export
DataAnalysis.__checkPointsPerYear <- function(startYear, endYear) {
  # check if for each point exactly one profile exists 
  # by checking if the "histograms" are equal (points per year == profiles per year)
  widths <- Loader.readInCurves(Paths.OSTALB_DATASET_WIDTHS, Symbols.REAL_DATA_SEPARATOR, TRUE);
  profiles <-  Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  
  perYearWidths <- CurvesMiner.collectWidthsData(widths)$widthsPerYear;
  perYearProfiles <- CurvesMiner.collectData(profiles)$profilesPerYear;

  # check
  print(Titles.DIFFERENCES_IN_YEARS);
  
  correct <- TRUE;
  for (i in startYear:endYear) {
    numberWidths <- length(perYearWidths[[toString(i)]]);
    numberProfiles <- length(perYearProfiles[[toString(i)]]);
    
    if (numberWidths != numberProfiles) {
      print(i);
      correct <- FALSE;
    }
  }
  
  if (correct) print(Strings.NONE);
}


#########################################################
#' Encodes data and sends it to a general histogram.
#'
#' @param startYear {number} the year at which the histogram should start
#' @param endYear {number} the year at which the histogram should end
#' @export
DataAnalysis.__createProfilesPerYear <- function(startYear, endYear) {
  years <- Encoder.createProfilesPerYearData(startYear, endYear);
  
  Plotter.plotGeneralHistogram(unlist(years), Titles.PROFILES_PER_YEAR, 
                               Strings.YEARS, Strings.FREQUENCY, TRUE, FALSE, 
                               Defaults.HISTOGRAM, FALSE);
}


#########################################################
#' Computes and stores the pairwise distances of the profiles from a year.
#'
#' @param startYear {number} the year at which to start
#' @param endYear {number} the year at which to end
#' @param normalize {logical} tells if the profiles should be normalized 
#' before computing a distance or not
DataAnalysis.__computePairwiseDistancesPerYear <- function(startYear, endYear, normalize = TRUE) {
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET, Symbols.REAL_DATA_SEPARATOR);
  perYearProfiles <- CurvesMiner.collectData(curvesData);
  
  profiles <- perYearProfiles$profilesPerYear;
  names <- perYearProfiles$namesPerYear;
  
  print(Titles.PROFILES_PER_YEAR);
  
  # iterate over all years
  for(i in startYear:endYear) {
    print(i);
    
    listOfProfiles <- profiles[[toString(i)]];
    listOfNames <- names[[toString(i)]];
    
    allDistancesPerProfile <- list();
    allNamesPerProfile <- list();
    
    # iterate over all profiles in that year
    for (j in 1:length(listOfProfiles)) {
      profile <- listOfProfiles[[j]];
      name <- listOfNames[[j]];
  
      distances <- Analyzer.__getPairwiseDistances(listOfProfiles[-(1:j)], profile, normalize);
      allDistancesPerProfile <- rlist::list.append(allDistancesPerProfile, distances);
      allNamesPerProfile <- rlist::list.append(allNamesPerProfile, name);
    }

    Storer.storeMatrix(allDistancesPerProfile, allNamesPerProfile, 
                       paste(i, Extensions.CSV, sep = Symbols.EMPTY));
  }
}


#########################################################
#' Creates an all-vs-all plot for the distances between profiles.
#'
#' @param startYear {number} the year at which to start
#' @param endYear {number} the year at which to end
#' @param meanDistances {logical} tells if the mean-distances should be plotted
#' @param title {string} the title for the plot
#' @export
DataAnalysis.__createAllVsAllDistancesPlot <- function(startYear, endYear, meanDistances = FALSE, 
                                                       title = Titles.PROFILES_DISTANCES_PER_YEAR) {
  distancesList <- Encoder.createAllVsAllDistancesPlotData(startYear, endYear, meanDistances, title);
  
  Plotter.plotGeneralBoxPlot(startYear:endYear, distancesList, title,
                             Strings.YEARS, Strings.DISTANCES, TRUE, FALSE);
}


#########################################################
#' Plots the ring-width consensus.
#' @param startYear {number} the year at which the check up should start
#' @param endYear {number} the year at which the check up should end
DataAnalysis.__createConsensus <- function(startYear, endYear) {
  # read in data
  curvesData <- Loader.readInCurves(Paths.OSTALB_DATASET_WIDTHS, Symbols.REAL_DATA_SEPARATOR, 
                                    TRUE, Extensions.META);

  # create consensus
  perYearWidths<- CurvesMiner.collectWidthsData(curvesData);
  masterChronology <- CurvesMiner.computeWidthsChronology(perYearWidths$widthsPerYear,
                                                          startYear,
                                                          endYear);
  
  # plot
  Plotter.plotGeneralLinePlot(unlist(masterChronology$years), unlist(masterChronology$widths), Titles.CONSENSUS, 
                              Strings.YEARS, Strings.WIDTH, TRUE, TRUE, xLimits = c(1910,2010));
}


#########################################################
# Moves through the chronology and compares neighboured consensus-profiles.
#'
#' @param consensiPath {string} the path to the consensi
#' @param consensusName {string} the name of the consensus
DataAnalysis.__compareCorrelationCoefficients  <- function(consensiPath, consensusName) {
  data <- Encoder.createCorrelationCoefficientHistogramEncoding(consensiPath, consensusName);
  
  xData <- data$xData;
  yData <- data$yData;
  
  Plotter.plotGeneralBarPlot(xData, yData, Titles.NEIGHBOURED_SCORES, 
                             Strings.YEARS, Strings.CORRELATION_COEFFICIENT, TRUE, FALSE,
                             yLimits = c(0.8, 1.0));
}