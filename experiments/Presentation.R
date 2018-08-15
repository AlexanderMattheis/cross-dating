#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

PRESENTATION_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Creates the visualized Presentation data.

#########################################################
#' Starts the ring width based cross-dating against density-based cross-dating.
Presentation.start <- function() {
  
  # Presentation.__createRankHistogram(Paths.PASSES_RING_WIDTH, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5), "p",
  #                                    Colors.GREEN, decreasing = TRUE);
  
  # Presentation.__createRankXyPlot(Paths.PASSES_RING_WIDTH, Files.SCORES_FOLDER_NAME,
  #                                 list(1, 2, 3, 4, 5), "p",
  #                                 Colors.GREEN, decreasing = TRUE);
  
  # Presentation.__createRankHistogram(Paths.PASSES, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5), "d",
  #                                    Colors.RED);
  
  # Presentation.__createRankXyPlot(Paths.PASSES, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5), "d",
  #                                    Colors.RED);
  
  # Presentation.__createRankXyPlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                 list(1, 2, 3, 4, 5), "d",
  #                                 Colors.BLUE);
  
  # Presentation.__createRankHistogram(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                    list(1, 2, 3, 4, 5), "d",
  #                                    Colors.BLUE);
  
  # Presentation.__createRankXyPlot(Paths.PASSES_RING_WIDTH_BUCKET, Files.SCORES_FOLDER_NAME,
  #                                 list(1, 2, 3, 4, 5), "p",
  #                                 Colors.VIOLET, decreasing = TRUE);
  
  # Presentation.__createPValueRankPlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                     list(1, 2, 3, 4, 5),  Defaults.TASK_D,
  #                                     Colors.BLUE, c(16, 15, 17, 18));
  
  # computing average runtime (ring-width approach)
  # Experiment4.__storeComputedScores(1, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  # Experiment4.__storeComputedScores(1, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  # Experiment4.__storeComputedScores(1, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  # Experiment4.__storeComputedScores(1, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  # Experiment4.__storeComputedScores(1, Defaults.TASK_PEARSON_CHAR, Defaults.CORRELATION_PEARSON);
  
  # computing average runtime (buckets approach)
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  # Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  
  # computing average runtime (buckets approach)
  # for (i in 1:5) {
  #   dates <- Interface.computeDatesHeuristicBucketApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #                                                          samplesPath = "input/interface/pass_1/samples/",
  #                                                          scoreTypeRingWidths = "p",
  #                                                          scoreTypeBuckets = "d",
  #                                                          topYearsCount = 20,
  #                                                          bestYearsMax = 5,
  #                                                          qualityMeasures = "s",
  #                                                          save = TRUE,
  #                                                          fileName = "datesBucketsHeuristic");
  # }
  
  # determine number of correct rated samples
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN_LENGTH_1, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5), "d",
  #                                         Palettes.NATURE);
  # 
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN_LENGTH_5, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5), "d",
  #                                         Palettes.NATURE);
  # 
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5), "d",
  #                                         Palettes.NATURE);
  # 
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_BUCKET_MIN_LENGTH_15, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5, 6, 7), "d",
  #                                         Palettes.NATURE);
  
  # draw mean rank curve
    # ring-width Pearson
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.099, 0.268, 0.459), Strings.QUALITY_2, Colors.GREEN);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.309, 0.552, 0.684), Strings.QUALITY_3, Colors.GREEN);
  
    # consensus method c)
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.049, 0.19, 0.316, 0.524), Strings.QUALITY_2, Colors.RED);
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.155, 0.406, 0.556, 0.688), Strings.QUALITY_3, Colors.RED);
  
    # per tree method c)
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.08, 0.364, 0.544, 0.688), Strings.QUALITY_2, Colors.DARK_BLUE);
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.231, 0.594, 0.76, 0.848), Strings.QUALITY_3, Colors.DARK_BLUE);
  
    # buckets method d)
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.081, 0.321, 0.588, 0.744), Strings.QUALITY_2, Colors.BLUE);
  # Presentation.plotWindowWidthDependency(c(1, 5, 10, 15), c(0.241, 0.606, 0.776, 0.879), Strings.QUALITY_3, Colors.BLUE);
  
    # voting method b)
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.305, 0.588, 0.74), Strings.QUALITY_2, Colors.DARK_RED);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.533, 0.752, 0.835), Strings.QUALITY_3, Colors.DARK_RED);
  
    # two-step Pearson + method d)
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.345, 0.58, 0.697), Strings.QUALITY_2, Colors.VIOLET);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.57, 0.728, 0.823), Strings.QUALITY_3, Colors.VIOLET);
  
  # create ring-width violine plot - compute ranks
  # Experiment1.createSummarizedViolinePlot(Paths.PASSES_RING_WIDTH_LENGTH_15, Files.SCORES_FOLDER_NAME,
  #                                         list(1, 2, 3, 4, 5, 6, 7),  list(Defaults.TASK_PEARSON_CHAR),
  #                                         c(Colors.ORANGE,
  #                                           Palettes.OTHER_NATURE[2],
  #                                           Palettes.OTHER_NATURE[3],
  #                                           Palettes.OTHER_NATURE[4]), TRUE);
  
  # measure runtimes
  # runtimes <- list();
  # for (i in 1:5) { # hint: deactivate console output before measuring
  #   # runtime <- Experiment4.__storeComputedScores(1, Defaults.TASK_T_VALUE_CHAR, Defaults.CORRELATION_PEARSON);
  #   # runtime <- Experiment1.__storeComputedScores(Defaults.TASK_C, 1, 0);
  #   # runtime <- Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, sum);
  #   # runtime <- Experiment2.storeComputedScores(Defaults.TASK_D, 1, min, identity);
  #   # runtime <- Experiment5.__computePerTreeScores(Defaults.TASK_C, 1);
  # 
  #   # measuring shifting runtime
  #   # runtime <- Interface.computeDatesVotingApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #   #                                                 samplesPath = "input/interface/pass_1/samples/",
  #   #                                                 scoreType = "d",
  #   #                                                 topYearsCount = 1,
  #   #                                                 approach = "p",  # activates powerset approach
  #   #                                                 minimumLength = 8,
  #   #                                                 bestYearsMax = 5,
  #   #                                                 save = TRUE,
  #   #                                                 fileName = "datesVoting");
  # 
  #   # runtime <- Interface.computeDatesHeuristicBucketApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #   #                                                        samplesPath = "input/interface/pass_1/samples/",
  #   #                                                        scoreTypeRingWidths = "p",
  #   #                                                        scoreTypeBuckets = "d",
  #   #                                                        topYearsCount = 20,
  #   #                                                        bestYearsMax = 5,
  #   #                                                        qualityMeasures = "s",
  #   #                                                        save = TRUE,
  #   #                                                        fileName = "datesBucketsHeuristic");
  # 
  #   # runtime <- Interface.computeDatesConsensusApproach(consensusPath = "input/interface/pass_1/",
  #   #                                                    consensusName = "consensus",
  #   #                                                    samplesPath = "input/interface/pass_1/samples/",
  #   #                                                    scoreType = "c",
  #   #                                                    bestYearsMax = 5,
  #   #                                                    save = TRUE,
  #   #                                                    fileName = "datesConsensus");
  # 
  #   # runtime <- Interface.computeDatesHeuristicBucketApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #   #                                                        samplesPath = "input/interface/pass_1/samples/",
  #   #                                                        scoreTypeRingWidths = "p",
  #   #                                                        scoreTypeBuckets = "d",
  #   #                                                        topYearsCount = 20,
  #   #                                                        bestYearsMax = 5,
  #   #                                                        qualityMeasures = "s",
  #   #                                                        save = TRUE,
  #   #                                                        fileName = "datesBucketsHeuristic");
  # 
  #   runtimes <- rlist::list.append(runtimes, runtime);
  # }
  # 
  #   print(unlist(runtimes));
  #   print(sum(unlist(runtimes))/5);
  
  # runtimes
  # Presentation.__createRuntimeDiagram(Approaches, c(0.5696873, 11.13604, 77.24634, 1721.60952, 200.61606, 279.0345));
  
    # two-step Pearson + method d) + maximum density
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.438, 0.696, 0.784), Strings.QUALITY_2, Colors.VIOLET);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.675, 0.836, 0.9), Strings.QUALITY_3, Colors.VIOLET);
  
    # maximum density Pearson
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.16, 0.448, 0.732), Strings.QUALITY_2, Colors.GREEN);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.505, 0.772, 0.952), Strings.QUALITY_3, Colors.GREEN);
  
    # ring-width t-value
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.066, 0.28, 0.554), Strings.QUALITY_2, Colors.GREEN);
  # Presentation.plotWindowWidthDependency(c(5, 10, 15), c(0.27, 0.504, 0.797), Strings.QUALITY_3, Colors.GREEN);
}


#########################################################
#' Creates a histogram out of the ranks from several computation passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} a scoreType
#' @param color {vector} the color you want use for the histogram
#' @param decreasing {logical} tells if the highest value should have lowest rank
#' @param exception {logical} tells if in the first round the decreasing parameter should be FALSE
#' @export
Presentation.__createRankHistogram <- function(passPathGenericName, scoresFolderGenericName, 
                                               passes, scoreType, color, decreasing = FALSE,
                                               exception = FALSE) {
  
  combinedRanksPerScoreType <- Encoder.createViolinePlotData(passPathGenericName, scoresFolderGenericName, 
                                                             passes, scoreType, decreasing,
                                                             exception); 
  yData <- combinedRanksPerScoreType[[1]];
  
  # (D): plot ranks
  Plotter.plotGeneralHistogram(yData[yData <= 20], Titles.RANKS, Strings.PREDICTED_CORRECT_YEAR_RANKS, Strings.FREQUENCY, 
                               TRUE, FALSE, Files.RANK_HISTOGRAM, TRUE, c(0, 160), color = color, xLimits = c(0,20));
}


#########################################################
#' Creates a xy-plot out of the ranks from several computation passes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} a scoreType
#' @param color {vector} the color you want use for the histogram
#' @param decreasing {logical} tells if the highest value should have lowest rank
#' @param exception {logical} tells if in the first round the decreasing parameter should be FALSE
#' @export
Presentation.__createRankXyPlot <- function(passPathGenericName, scoresFolderGenericName, 
                                            passes, scoreType, color, decreasing = FALSE,
                                            exception = FALSE) {
  
  combinedRanksPerScoreType <- Encoder.createViolinePlotData(passPathGenericName, scoresFolderGenericName, 
                                                             passes, scoreType, decreasing,
                                                             exception); 
  xData <- 1:20;
  yData <- c(combinedRanksPerScoreType[[1]], 1:20);  # or zero bars not added for the first 20
  yData <- as.vector(table(yData));
  yData <- yData - 1;  # above it was added one elemnt to the first 20 bars, now it is removed
  
  # (D): plot ranks
  Plotter.plotGeneralLinePlot(xData = xData, 
                              yData =  yData[1:20], 
                              plotTitle = Titles.RANKS, 
                              xAxisLabel = Strings.PREDICTED_CORRECT_YEAR_RANKS, 
                              yAxisLabel = Strings.FREQUENCY, 
                              print = TRUE, 
                              save = FALSE,
                              points = TRUE, 
                              xLimits = c(0,20), 
                              yLimits = c(0, 80),
                              fitCurve = FALSE, 
                              percentage = FALSE, 
                              color = color);
}


#########################################################
#' Creates a line plot showing the dependency 
#' between length and mean ranks for example.
#'
#' @param xData {vector} the data for the x-axis
#' @param yData {vector} the corresponding data for the y-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param color {string} the color used for the curve
#' @export
Presentation.plotWindowWidthDependency <- function(xData, yData, yAxisLabel, color) {
  Plotter.plotGeneralLinePlot(xData, yData, plotTitle = Titles.WIDNOW_QUALITY, 
                              xAxisLabel = Strings.WINDOW_LENGTH, yAxisLabel = yAxisLabel, TRUE, FALSE,
                              points = TRUE, xLimits = c(0, 15), yLimits = c(0, 1.0), percentage = TRUE, color = color);
}


#########################################################
#' Creates p-value/rank plot. So an xy plot with the given axes.
#'
#' @param passPathGenericName {string} the generic name of the pass path (without the pass number)
#' @param scoresFolderGenericName {string} the generic name of the score type path within pass path (without the scoreType)
#' @param passes {list} contains the numbers from different passes which should be used for computation
#' @param scoreType {string} the scoreType for which the plot should be created
#' @param palette {vector} the colors palette which should be used
#' @param shapes {vector} the shapes palette which should be used
#' @param weighted {logical} tells if the lowest score should be weighted with its delta score or not 
#'
#' @example 
#' x = delta score = PredictionRank2Score-PredictionRank1Score (because PredictionRank2Score > PredictionRank1Score), 
#' y = (real) rank
Presentation.__createPValueRankPlot <- function(passPathGenericName, scoresFolderGenericName, 
                                                passes, scoreType, palette, shapes, weighted = FALSE) {
  
  data <- Encoder.createRankPValueEncoding(passPathGenericName, scoresFolderGenericName, 
                                           passes, scoreType, 1);
  
  xData <- data$pValuesPerSample;
  yData <- data$ranksPerSample;
  group <- rep(Defaults.P_VALUE_GROUP, length(yData));
  
  # print(Titles.P_VALUES_BELOW_1_PERCENT);
  # print(length(xData[xData <= 0.01]));
  # 
  # print(Titles.P_VALUES_BELOW_0_1_PERCENT);
  # print(length(xData[xData <= 0.001]));
  
  orderOfScores <- order(xData);
  sortedRanks <- yData[orderOfScores];
  firstHundredRanks <- sortedRanks[1:100];
  
  print(Strings.NUMBER_RANK_1);
  print(length(firstHundredRanks[firstHundredRanks == 1]));
  
  Plotter.plotGeneralScatterPlot(unlist(xData), unlist(yData), unlist(group), 
                                 Titles.RANKS_FOR_P_VALUES, 
                                 Defaults.P_VALUES,
                                 Strings.PREDICTED_CORRECT_YEAR_RANKS, Strings.RANKS, palette, shapes,
                                 TRUE, FALSE, xLimits = c(0, 0.105),
                                 xLogarithmic = TRUE);
}


#########################################################
#' Creates a runtime diagram.
Presentation.__createRuntimeDiagram <- function(names, runtimes) {
  Plotter.plotGeneralBarPlot(xData = names, yData = runtimes, plotTitle = Strings.RUNTIMES, 
                             xAxisLabel = Strings.APPROACH, yAxisLabel = Strings.RUNTIME, 
                             print = TRUE, save = FALSE, flip = TRUE);
}