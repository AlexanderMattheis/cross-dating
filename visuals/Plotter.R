#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

PLOTTER_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Contains several methods to visualize the data.

#########################################################
#' Plots the distribution of given scores in a histogram.
#' On the x-axis you see the different scores and on the y-axis the frequency in percent.
#'
#' @param scoreListsPerPatternLength {list} output of the "computeMicaScores"-function
#' @param correctPosition {numerical} the correct position for the pattern which was aligned against every position
#' @param save {logical} save the plot as pdf or not
#' @param print {logical} show the plot or not
#' @param parametersToEncode {list} the list of parameters which should be encoded in the filename
#'
#' @export
Plotter.plotHistogramsForScores <- function(scoreListsPerPatternLength, correctPosition, save, print, parametersToEncode) {
  xAxisMarks <- Encoder.getCorrectPatternPositionScores(scoreListsPerPatternLength, correctPosition, FALSE);  # scores of correct position
  
  # iterate over all pattern lenghts
  for (i in 1:length(scoreListsPerPatternLength)) {
    scores <- scoreListsPerPatternLength[[i]];
    values <- as.numeric(scores);
    
    # hint: "y" has to be changed, if you change "binwidth"
    histogram <- qplot(x = values, y = 0.1*..density.., xlab = Strings.SCORE, ylab = Strings.FREQUENCY, 
                       main = paste(paste(Titles.SCORE_FREQUENCY_DIAGRAM_START, toString(i), sep = Symbols.EMPTY), 
                                    Titles.SCORE_FREQUENCY_DIAGRAM_END, sep = Symbols.EMPTY),
                       geom = "histogram", fill = I(Colors.GREEN), col = I(Colors.WHITE), binwidth = 0.1);
    
    plot <- histogram + 
      scale_y_continuous(labels = scales::percent) +
      geom_vline(xintercept = xAxisMarks[[i]], linetype = "dashed", color = Colors.DARK_BROWN, size = 1) + 
      geom_density(alpha = 0.2, fill = Colors.RED, color = alpha(Colors.RED, 0.1), size = 0);
    
    if (print) {
      print(plot);
    }
    
    if (save) {
      Storer.saveEncodedPlot(plot, i, parametersToEncode);
    }
  }
}


#########################################################
#' Creates for each year in first curve and second curve a common plot for the two profiles.
#'
#' @param curve1Data {list(years, profiles, name)} the unique years and the profiles
#' @param curve2Data {list(years, profiles, name)} the unique years and the profiles
#' @param plotTitle {string} the title of the plot
#' @param legendTitle {string} the title of the legend
#' @param curve1Title {string} the title of the first curve
#' @param curve2Title {string} the title of the second curve
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as pdf or not
#' @param fileExtension {string} the type of the stored file
#' @param profilesNumPerYear {list(masterChronology, consensus)} stores the number of profiles per Year
#' @param yLimits {vector} limits for the y-axis
Plotter.plotProfiles <- function(curve1Data, curve2Data, plotTitle, legendTitle, 
                                 curve1Title, curve2Title, print, save, fileExtension, 
                                 profilesNumPerYear, yLimits = vector()) {
  # check if both curves have same years and same number of years
  years1 <- curve1Data$years;
  years2 <- curve2Data$years;
  
  if (!identical(years1, years2)) stop(Exceptions.CURVE_YEARS_NOT_IDENTICAL);
  if (length(profilesNumPerYear$masterChronology) != 
      length(profilesNumPerYear$consensus)) stop(Exceptions.CURVES_NUM_PROFILES);
  
  # iterate over each year
  for (i in 1:length(years1)) {
    # encode data for ggplot
    profile1 <- curve1Data$profiles[[i]];
    profile2 <- curve2Data$profiles[[i]];
    
    data1 <- data.frame(x1 = (1:length(profile1)) / length(profile1), y1 = unlist(profile1));
    data2 <- data.frame(x2 = (1:length(profile2)) / length(profile2), y2 = unlist(profile2));
    # plot
    plot <- ggplot() + 
      geom_line(data = data1, aes(x = x1, y = y1, color = Colors.ORANGE), size = 2.5) + 
      geom_line(data = data2, aes(x = x2, y = y2, color = Colors.GREEN), size = 1.5) +
      labs(
        title = paste(plotTitle, years1[[i]]),
        x = Strings.TIME,
        y = Strings.DENSITY
      ) +
      scale_colour_manual(name = legendTitle, 
                          labels = c(paste(curve1Title, Symbols.SPACE, Symbols.BRACKET_LEFT, 
                                           profilesNumPerYear$consensus[[i]], 
                                           Symbols.BRACKET_RIGHT, sep = Symbols.EMPTY),
                                     paste(curve2Title, Symbols.SPACE, Symbols.BRACKET_LEFT, 
                                           profilesNumPerYear$masterChronology[[i]], 
                                           Symbols.BRACKET_RIGHT, sep = Symbols.EMPTY)),
                          values = c(Colors.ORANGE, Colors.GREEN));
    
    if (length(yLimits) > 0) {
      plot <- plot + coord_cartesian(ylim = yLimits); 
    }
    
    if (print) print(plot);
    if (save) Storer.savePlot(plot, paste(Defaults.CONSENSUS_PLOT_PDF_NAME, 
                                          years1[[i]], fileExtension, sep = Symbols.EMPTY));
  }
}


#########################################################
#' Plots the rank-violine plots of different score-types side by side.
#'
#' @param ranksPerScoreType {list} list of sublists, 
#' where every sublist contains the ranks of the different samples
#' @param types {list} list of possible score types 
#' (maximum number of different score types = length(palette) and ggplot internal maximum)
#' @param plotTitle {string} the title of the plot
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param fileExtension {string} the type of the stored file
#' @param palette {vector} the colors palette which should be used
#' @param xLabel {string} the label for the x-axis
#' @param printNumRanks {logical} tells if it should be printed out the number of with rank x rated samples
#' 
#' @return {list(median, mean, variance)} The computed values for the last distribution
#' @export
Plotter.plotViolinePlots <- function(ranksPerScoreType, types, plotTitle, 
                                     print, save, fileExtension, palette, 
                                     xLabel = Strings.TYPE, printNumRanks = TRUE) {
  plot <- ggplot();
  
  yMedian <- 0;
  yVariance <- 0;
  yMean <- 0;
  
  if (length(types) <= length(palette)) {  # check if enough colors to plot
    # iterate over each scoreType
    for (i in 1:length(ranksPerScoreType)) {
      # encode data
      xType <- types[[i]];
      yRanks <- unlist(ranksPerScoreType[[i]]);
      
      yMedian <- median(yRanks);
      yVariance <- round(var(yRanks), 2);
      yMean <- round(mean(yRanks), 2);

      if (printNumRanks) {
        print(xType);
        Plotter.__printNumRanks(yRanks, Parameters.RANKS_TO_OUTPUT);
      }
      
      # tests
      # if (is.na(yMedian)) {
      #   print(i);
      # }
      
      xInformation <- paste(xType, 
                            paste(Strings.MEDIAN, yMedian), 
                            paste(Strings.MEAN, yMean), 
                            paste(Strings.VARIANCE, yVariance), 
                            sep = Symbols.NEW_LINE);
      
      data <- data.frame(type = xInformation, ranks = yRanks);
      
      plot <- plot + 
        geom_violin(data = data, aes(x = type, y = ranks), 
                    fill = palette[[i]],
                    color = Colors.WHITE, size = 1) +
        labs(
          title = plotTitle,
          x = xLabel,
          y = Strings.RANK
        ) + 
        theme(axis.text.x = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0,80));
    }
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.VIOLINE_PLOT_NAME, 
                                           fileExtension, sep = Symbols.EMPTY));
  
  return(list(median = yMedian, mean = yMean, variance = yVariance));
}


#########################################################
#' Prints number of ranks with given number of ranks. 
#' Minimally the number of with rank 1 rated samples is printed out.
#'
#' @param yRanks {vector} the ranks in which it is counted
#' @param maxRankToCheck {number} tells the maximum rank for which the number should be printed
Plotter.__printNumRanks <- function(yRanks, maxRankToCheck) {
  allRated <- length(yRanks);
  
  for (i in 1:maxRankToCheck) {
    numRatedSamples <- length(which(yRanks == i));
    print(paste(Strings.NUM_RANK_, i, Symbols.SPACE, 
                numRatedSamples, Symbols.SLASH, allRated, sep = Symbols.EMPTY));
  }

  rankBiggerRated <- length(which(yRanks > maxRankToCheck));
  
  print(paste(Strings.NUM_RANK_BIGGER, Symbols.SPACE, 
              rankBiggerRated, Symbols.SLASH, allRated, sep = Symbols.EMPTY));
}


#########################################################
#' Creates a side by side plot for testing purposes 
#' of the curves to align and their consensus.
#'
#' @param profilesToAlign {dataframe} the curves which are aligned together
#' @param consensusCurve {list} x and y coordinates from a consensus-curve
#' @param year {numerical} current year
#' @param yLimits {vector} limits for the y-axis
Plotter.plotConsensusAndCurves <- function(profilesToAlign, consensusCurve, year, yLimits = vector()) {
  # curves plot
  profilesToAlign$x <- (1:nrow(profilesToAlign)) / nrow(profilesToAlign);
  plotData <- reshape2::melt(profilesToAlign, id.var = "x");
  curvesPlot <- ggplot(plotData, aes(x = x, y = value, group = variable, colour = variable)) +
    labs(
      title = Titles.CURVES,
      x = Strings.TIME,
      y = Strings.DENSITY
    ) +
    theme(legend.position = "none") +
    geom_line(aes(lty=variable), size = 0.5);
  
  # consensus plot
  consensus <- data.frame(consensusX = consensusCurve$x, consensusY = consensusCurve$y);
  consensusPlot <- ggplot(data = consensus, aes(x = consensusX, y = consensusY)) +
    labs(
      title = Titles.CONSENSUS,
      x = Strings.TIME,
      y = Strings.DENSITY
    ) +
    geom_line(color = Colors.GREEN, size = 0.5);

  if (length(yLimits) > 0) {
    curvesPlot <- curvesPlot + coord_cartesian(ylim = yLimits); 
    consensusPlot <- consensusPlot + coord_cartesian(ylim = yLimits);
  }
  
  # combined plot
  combinedPlot <- gridExtra::grid.arrange(curvesPlot, consensusPlot, ncol=2);
  
  Storer.savePlot(combinedPlot, paste(toString(year), Extensions.PDF, sep = Symbols.EMPTY));
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
#' @export
Plotter.plotVennDiagram <- function(ranksPerScoreType, fileNames, ranksToLookAt, pass, save) {
  fileNamesPerSet <- Encoder.creatVennDiagram1Data(ranksPerScoreType, fileNames, ranksToLookAt, pass, save);
  
  Plotter._drawVennDiagram(fileNamesPerSet, pass, save);
}


#########################################################
#' Draw the Venn Diagram given a set of different sets to display.
#'
#' @param dataPerSet {list} the sets to display
#' @param pass {numerical} the pass to encode into file names 
#' @param save {logical} save plot and a list with maximum intersection elements or not
Plotter._drawVennDiagram <- function(dataPerSet, pass, save) {
  # retrieve information for plot #
  # compute set sizes
  area1 <- length(dataPerSet[[1]]);
  area2 <- length(dataPerSet[[2]]);
  area3 <- length(dataPerSet[[3]]);
  area4 <- length(dataPerSet[[4]]);
  
  # compute bi-intersections
  intersection12 <- Reduce(intersect, dataPerSet[c(1,2)]);
  intersection13 <- Reduce(intersect, dataPerSet[c(1,3)]);
  intersection14 <- Reduce(intersect, dataPerSet[c(1,4)]);
  
  intersection23 <- Reduce(intersect, dataPerSet[c(2,3)]);
  intersection24 <- Reduce(intersect, dataPerSet[c(2,4)]);
  
  intersection34 <- Reduce(intersect, dataPerSet[c(3,4)]);
  
  # compute tri-intersections
  intersection123 <- Reduce(intersect, dataPerSet[c(1,2,3)]);
  intersection124 <- Reduce(intersect, dataPerSet[c(1,2,4)]);
  intersection134 <- Reduce(intersect, dataPerSet[c(1,3,4)]);
  
  intersection234 <- Reduce(intersect, dataPerSet[c(2,3,4)]);
  
  # compute quatro-intersections
  intersection1234 <- Reduce(intersect, dataPerSet[1:4]);
  
  # draw
  venn.plot <- draw.quad.venn(area1 = area1,
                              area2 = area2,
                              area3 = area3,
                              area4 = area4,
                              n12 = length(intersection12),
                              n13 = length(intersection13),
                              n14 = length(intersection14),
                              n23 = length(intersection23),
                              n24 = length(intersection24),
                              n34 = length(intersection34),
                              n123 = length(intersection123),
                              n124 = length(intersection124),
                              n134 = length(intersection134),
                              n234 = length(intersection234),
                              n1234 = length(intersection1234),
                              alpha = c(0.5, 0.5, 0.5, 0.5),
                              category = Parameters.SCORE_TYPES_TO_LOOK_AT,
                              cat.fontfamily = rep("Helvetica", 4),
                              col = Palettes.NATURE[1:4],
                              cat.cex = rep(1.25, 4),
                              fill = Palettes.NATURE[1:4],
                              fontfamily = rep("Helvetica", 15),
                              label.col = Colors.WHITE,
                              lwd = "blank");
  
  if (save) { 
    Storer.saveVennDiagram(intersection1234, pass); 
  }
}


#########################################################
#' Creates a Venn diagram for the ranks of 4 measurement types.
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
#' @export
Plotter.plotVennDiagram2 <- function(ranksPerScoreType, fileNames, rankIndicesToLookAt, pass, save) {
  fileNamesPerSet <- Encoder.creatVennDiagram2Data(ranksPerScoreType, fileNames, rankIndicesToLookAt, pass, save);
  
  # draw
  Plotter._drawVennDiagram(fileNamesPerSet, pass, save);
}


#########################################################
#' For different score-thresholds d in [d_min, d_max]
#' it is looked how many ranks are in the top t ranks
#' or not there under the given threshold d i.e. 
#' every score above the thresold is cut away.
#' The two numbers of how many there (recall) or not there (fallout)
#' are then plotted in percentage into the same plot,
#' where as for the recall the y-axis is used
#' and for the fallout the x-axis.
#' Hint: If there is only a single scoreType, then d-values are plotted.
#'
#' @param dataPerScoreType {list(ranksPerScoreType, scoresPerScoreType)} two lists of sublists, 
#' where every sublist contains the ranks (respectively the score) of the different samples
#' @param types {list} list of possible score types 
#' (maximum number of different score types = length(palette) and ggplot internal maximum)
#' @param plotTitle {string} the title of the plot
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param fileExtension {string} the type of the stored file
#' @param maxTopRank {numerical} the maximum number of top interval [1, 2, .., maxTopRank]
#' @param color {string} the color if there is only a single type
#' 
#' @return {list(falsePositives, truePositives)} the list of false and true positives for unit-testing
#' @export
Plotter.plotROC <- function(dataPerScoreType, types, plotTitle, print, save, 
                            fileExtension, maxTopRank = Parameters.TOP_RANK,
                            color = Colors.BLUE) {

  data <- Encoder.createRocPlotEncoding(dataPerScoreType, types, maxTopRank);
  
  falsePositives <- data$falsePositives;
  truePositives <- data$truePositives;
  tresholdValues <- data$tresholdValues;
  type <- data$type;
  
  # check if only one type
  plotScores <- length(types) == 1;
  
  # conversion
  types <- paste(types, Symbols.BRACKET_RIGHT, sep = Symbols.EMPTY);
  
  # plot
  data <- data.frame(g = factor(unlist(type)), x = unlist(falsePositives), y = unlist(truePositives));
  
  plot <- ggplot2::ggplot(data = data, aes(x = x, y = y)) + 
    geom_line(aes(color = g)) + 
    labs(
      title = plotTitle,
      x = Strings.FALSE_POSITIVE_RATE,
      y = Strings.TRUE_POSITIVE_RATE
    ) + geom_abline(color = Colors.GRAY, intercept = 0, slope = 1);
  
  if (plotScores) {  # add scores to graph
    n <- round(length(falsePositives) / length(types));
    n12 <- round(n/2);
    n14 <- round(n/4);
    n34 <- round((3*n)/4);
    
    plot <- plot + 
      scale_color_manual(name = Titles.METHOD, labels = types, values = color) + 
      geom_point(aes(x = x[1], y = y[1]), color = color) +
      geom_point(aes(x = x[n14], y = y[n14]), color = color) +
      geom_point(aes(x = x[n12], y = y[n12]), color = color) +
      geom_point(aes(x = x[n34], y = y[n34]), color = color) +
      geom_point(aes(x = x[n], y = y[n]), color = color) +
      geom_text(aes(x = x[1], y = y[1], label = unlist(tresholdValues[1])), color = color, nudge_x = 0.07, nudge_y = 0.03) +
      geom_text(aes(x = x[n14], y = y[n14], label = unlist(tresholdValues[n14])), color = color, nudge_x = 0.07, nudge_y = -0.04) +
      geom_text(aes(x = x[n12], y = y[n12], label = unlist(tresholdValues[n12])), color = color, nudge_x = 0.07, nudge_y = -0.04) +
      geom_text(aes(x = x[n34], y = y[n34], label = unlist(tresholdValues[n34])), color = color, nudge_x = 0.07, nudge_y = -0.04) +
      geom_text(aes(x = x[n], y = y[n], label = unlist(tresholdValues[n])), color = color, nudge_x = -0.04, nudge_y = -0.04);
  } else {
    plot <- plot + scale_color_manual(name = Titles.METHOD, labels = types, values = Palettes.NATURE[1:length(types)]);
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.ROC_PLOT_NAME, 
                                        fileExtension, sep = Symbols.EMPTY));
  
  return(list(truePositives = truePositives, falsePositives = falsePositives));
}


#########################################################
#' Creates a back-to-back histogram from given ranks.
#'
#' @param ranksPerScorePass {list} list of ranks per score type
#' @param fileNamesPerPass {list} the names fo the files per pass
#' @param rankIndicesToLookAt {vector} the rank indices in ascending ordered ranks for which the Venn Diagrams should be created 
#' (negative ranks indices gives you worst ranks i.e. -1 = last rank)
#' @param maxIntersection {logical} if you want create a plot for the intersection of all four methods
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @export
Plotter.plotBackToBackHistogram <- function(ranksPerScorePass, fileNamesPerPass, 
                                            rankIndicesToLookAt, maxIntersection, print, save) {
  # retrieve data
  positiveYears <- Encoder.getYearHistogramEncoding(ranksPerScorePass, fileNamesPerPass,
                                                    rankIndicesToLookAt, maxIntersection);
  
  negativeYears <- Encoder.getYearHistogramEncoding(ranksPerScorePass, fileNamesPerPass, 
                                                    -rankIndicesToLookAt, maxIntersection);
  
  # encode
  length(positiveYears) <- max(length(positiveYears), length(negativeYears));
  length(negativeYears) <- max(length(positiveYears), length(negativeYears));
  data = data.frame(positive = positiveYears, negative = negativeYears);
  
  # create
  histogramPositive <- geom_histogram(data = data, mapping = aes(x = positive, y = ..count..), 
                                      fill = I(Colors.GREEN), col = I(Colors.WHITE), binwidth = 1);

  histogramNegative <- geom_histogram(data = data, mapping = aes(x = negative, y = -..count..),
                                      fill = I(Colors.RED), col = I(Colors.WHITE), binwidth = 1);
  
  # plot
  plot <- ggplot(data, aes(x = positive)) + 
    histogramNegative + histogramPositive + 
    labs(
      title = Titles.DISTRIBUTION_OF_YEARS,
      x = Strings.YEARS,
      y = Strings.FREQUENCY
    );
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.BACK_TO_BACK_HISTOGRAM, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
}


#########################################################
#' Creates a histogram with the given data.
#'
#' @param data {vector} the vector containing the elements which have to be counted
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param fileName {string} the name for the file
#' @param barLabels {logical} tells if to plot lables above bars or not
#' @param yLimits {vector} limits for the y-axis (necessary if barLables)
#' @param xLimits {vector} limits for the x-axis
#' @param binwidth {numeric} the width of the bins
#' @param drawDensity {logical} tells if a density-curve should be drawn
#' @param percentage {logical} tells if teh y-axis should use percent or not
#' @param color {string} the color you want use
#' 
#' @return {ggplot} the plot
#' @export
Plotter.plotGeneralHistogram <- function(data, plotTitle, xAxisLabel, yAxisLabel, 
                                         print, save, fileName, barLabels = TRUE, yLimits = vector(), xLimits = vector(),
                                         binWidth = 1, drawDensity = FALSE, percentage = FALSE, color = Colors.GREEN) {
  # encode
  data = data.frame(x = data);
  
  # create
  plot <- ggplot(data, aes(x = x)) +
    geom_histogram(data = data, mapping = aes(x = x, y = ..count..), 
                   fill = I(color), col = I(Colors.WHITE), binwidth = binWidth) +
    labs(
      title = plotTitle,
      x = xAxisLabel,
      y = yAxisLabel
    );
  
  if (barLabels) {
    plot <- plot + stat_bin(binwidth= binWidth, geom="text", aes(label=..count..), vjust = -1);
  }
  
  if (length(yLimits) != 0) {  # if vector not empty
    plot <- plot + coord_cartesian(ylim = yLimits);
  }
  
  if (length(xLimits) != 0) {
    plot <- plot + coord_cartesian(xlim = xLimits);
  }
  
  if (length(xLimits) != 0 && length(yLimits) != 0) {
    plot <- plot + coord_cartesian(xlim = xLimits, ylim = yLimits);
  }
  
  if (drawDensity) {
    plot <- plot + geom_density(aes(y = ..count.. * 1));
  }
  
  if (percentage) {
    plot <- plot + scale_y_continuous(labels = scales::percent);
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(fileName, Extensions.PNG, sep = Symbols.EMPTY));
  
  return(plot);
}


#########################################################
#' Creates a box plot with the given data.
#'
#' @param xData {vector} the vector with x-axis data
#' @param yData {list} the list of vectors at a specific x-coordinate
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @export
Plotter.plotGeneralBoxPlot <- function(xData, yData, plotTitle, xAxisLabel, yAxisLabel, print, save) {
  # encode data
  plot <- ggplot();
  
  for (i in 1:length(xData)) {
      # encode data
      xType <- rep(xData[[i]], length(yData[[i]]));
      yType <- unlist(yData[[i]]);

      data <- data.frame(x = xType, y = yType);
      
      plot <- plot + 
        geom_boxplot(data = data, aes(x = x, y = y), 
                     outlier.colour = Colors.RED,
                     fill = Colors.GREEN,
                     color = Colors.WHITE, size = 0.5) +
        labs(
          title = plotTitle,
          x = xAxisLabel,
          y = yAxisLabel
        );
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.BOX_PLOT_NAME, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
}


#########################################################
#' Creates a bar plot with the given data.
#'
#' @param xData {vector} the vector with x-axis data
#' @param yData {vector} the vector with y-axis data
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param barLabels {logical} tells if counts should be put above the bars
#' @param yLimits {vector} limits for the y-axis (necessary if barLables)
#' @param centerTitle {logical} tells if the title should be placed to the center
#' @param color {string} the string for the fill of a bar
#' @param flip {logical} tells if the bar should be drawn horizontly
#' 
#' @return {ggplot} the plot
#' @export
Plotter.plotGeneralBarPlot <- function(xData, yData, plotTitle, xAxisLabel, yAxisLabel, print, save, 
                                       barLabels = FALSE, yLimits = vector(), centerTitle = FALSE, 
                                       color = Colors.GREEN, flip = FALSE) {
  # encode data
  data <- data.frame(x = xData, y = yData);
  
  # create plot
  plot <- ggplot() + geom_bar(stat="identity", data = data, 
                              aes(x = x, y = y), 
                              fill = color,
                              color = Colors.WHITE, 
                              size = 0.5) +
    labs(
      title = plotTitle,
      x = xAxisLabel,
      y = yAxisLabel
    );  
  
  # additions
  if (barLabels) {  
    plot <- plot + geom_text(data = data, aes(x = x, y = y, label = y), vjust = -0.5);
  }
  
  if (length(yLimits) != 0) { # if vector not empty
    plot <- plot + coord_cartesian(ylim = yLimits);
  }
  
  if (centerTitle) {
    plot <- plot + theme(plot.title = element_text(hjust = 0.5));
  }
  
  if (flip) {
    plot <- plot + coord_flip() + scale_x_discrete(limits = xData);  # to avoid resorting
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.BAR_PLOT_NAME, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
  
  return(plot);
}


#########################################################
#' Plots the rank-violine plots of different score-types side by side.
#'
#' @param ranksPerScoreType {list} list of sublists, 
#' where every sublist contains the ranks of the different samples
#' @param types {list} list of possible score types 
#' (maximum number of different score types = length(palette) and ggplot internal maximum)
#' @param plotTitle {string} the title of the plot
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param fileExtension {string} the type of the stored file
#' @param palette {vector} the colors palette which should be used
#' @param xLabel {string} the label for the x-axis
#' @param printNumRanks {logical} tells if it should be printed out the number of with rank x rated samples
#' 
#' @return {list(median, mean, variance)} The computed values for the last distribution
#' @export
Plotter.plotBoxPlots <- function(ranksPerScoreType, types, plotTitle, 
                                 print, save, fileExtension, palette, 
                                 xLabel = Strings.TYPE, printNumRanks = TRUE) {
  plot <- ggplot();
  
  yMedian <- 0;
  yVariance <- 0;
  yMean <- 0;
  
  if (length(types) <= length(palette)) {  # check if enough colors to plot
    # iterate over each scoreType
    for (i in 1:length(ranksPerScoreType)) {
      # encode data
      xType <- types[[i]];
      ranks <- unlist(ranksPerScoreType[[i]]);
      yRanks <- ranks[ranks > 0];
      outliers <- ranks[ranks < 1];

      yMedian <- median(yRanks);
      yVariance <- round(var(yRanks), 2);
      yMean <- round(mean(yRanks), 2);
      
      if (printNumRanks) {
        print(xType);
        Plotter.__printNumRanks(yRanks, 5);
      }
      
      # tests
      # if (is.na(yMedian)) {
      #   print(i);
      # }
      
      xInformation <- paste(xType,
                            paste(Strings.OUTLIERS, length(outliers)),
                            paste(Strings.MEDIAN, yMedian), 
                            paste(Strings.MEAN, yMean), 
                            paste(Strings.VARIANCE, yVariance), 
                            sep = Symbols.NEW_LINE);

      data <- data.frame(type = xInformation, ranks = yRanks);
      
      plot <- plot + 
        geom_boxplot(data = data, aes(x = type, y = ranks), 
                     fill = palette[[i]],
                     color = Colors.WHITE, size = 0.2) +
        labs(
          title = plotTitle,
          x = xLabel,
          y = Strings.PEAK_RANK
        ) + 
        theme(axis.text.x = element_text(hjust = 0.5));
    }
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.BOX_PLOT_NAME, 
                                        fileExtension, sep = Symbols.EMPTY));
  
  return(list(median = yMedian, mean = yMean, variance = yVariance));
}


#########################################################
#' Plots a custom dot plot.
#' 
#' @param xData {vector} the vector with x-axis data
#' @param yData {vector} the vector with y-axis data
#' @param group {vector} the vector with the group number of the y-value
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param groupLabel {string} the label for the group
#' @param colors {vector} the colors to use
#' @param shapes {vector} the shapes to use
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param sizes {vector} changes the sizes of plotted elements
#' @param xLimits {vector} the 2D vector specifying upper and lower bound
#' @param yLimits {vector} the 2D vector specifying upper and lower bound
#' @param xLogarithmic {logical} tells if the x-axis should be logarithmic or not
#' 
#' @return {ggplot} the plot
#' @export
Plotter.plotGeneralScatterPlot <- function(xData, yData, group, plotTitle, 
                                           xAxisLabel, yAxisLabel, groupLabel, colors, shapes,
                                           print, save, sizes = NULL, 
                                           xLimits = vector(), yLimits = vector(), xLogarithmic = FALSE) {
  # encode data
  data <- data.frame(x = xData, y = yData, g = group);
  
  # create plot
  if (!is.null(sizes)) {
    plot <- ggplot(data, aes(x = x, y = y, shape = g, color = g, size = g)) + 
      scale_size_manual(values = sizes) + 
      labs(
        title = plotTitle,
        x = xAxisLabel,
        y = yAxisLabel,
        color = groupLabel,
        shape = groupLabel,
        size = groupLabel
      );
  } else {
    plot <- ggplot(data, aes(x = x, y = y, shape = g, color = g)) + 
      labs(
        title = plotTitle,
        x = xAxisLabel,
        y = yAxisLabel,
        color = groupLabel,
        shape = groupLabel
      );
  }
  
  plot <- plot + geom_point() +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    theme(legend.position = "none");
  
  # set limits
  if (length(xLimits) > 0 && length(yLimits) > 0) {
    plot <- plot + coord_cartesian(xlim = xLimits, ylim = yLimits);
  } else if (length(xLimits) > 0) {
    plot <- plot + coord_cartesian(xlim = xLimits);
  } else if (length(yLimits) > 0) {
    plot <- plot + coord_cartesian(ylim = yLimits);
  }
  
  if (xLogarithmic) {
    plot <- plot + coord_trans(x = "log10") + scale_x_continuous(breaks = c(0.00001, 0.001, 0.01, 0.1));
  }
  
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.SCATTER_PLOT_NAME, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
  
  return(plot);
}


#########################################################
#' Plots a custom line plot.
#'
#' @param xData {vector} the vector with x-axis data
#' @param yData {vector} the vector with y-axis data
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param points {logical} tells if points should be plotted or not
#' @param xLimits {vector} the 2D vector specifying upper and lower bound
#' @param yLimits {vector} the 2D vector specifying upper and lower bound
#' @param fitCurve {logical} tells if a curve should be fitted through points or not
#' @param percentage {logical} tells if on the y-axis are percentages or not 
#' @param xBreaks {numeric} tells how many breaks should be created on the x-axis
#' 
#' @return {ggplot} the plot
#' @export
Plotter.plotGeneralLinePlot <- function(xData, yData, plotTitle, 
                                        xAxisLabel, yAxisLabel, print, save,
                                        points = FALSE, xLimits = vector(), yLimits = vector(),
                                        fitCurve = FALSE, percentage = FALSE, color = Colors.GREEN, xBreaks = -1) {
  # encode data
  data <- data.frame(x = xData, y = yData);
  
  # create plot
  plot <- ggplot(data, aes(x = x, y = y)) + geom_line(color = color, size = 0.5) +
    labs(
      title = plotTitle,
      x = xAxisLabel,
      y = yAxisLabel
    );
  
  # add points
  if (points) {
    plot <- plot + geom_point(color = color);
  }
  
  if (fitCurve) {
    plot <- plot + geom_smooth();
  }
  
  # set limits
  if (length(xLimits) > 0 && length(yLimits) > 0) {
    plot <- plot + coord_cartesian(xlim = xLimits, ylim = yLimits);
  } else if (length(xLimits) > 0) {
    plot <- plot + coord_cartesian(xlim = xLimits);
  } else if (length(yLimits) > 0) {
    plot <- plot + coord_cartesian(ylim = yLimits);
  }
  
  if (percentage) {
    plot <- plot + scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 6));
  }
  
  if (xBreaks >= 0) {
    plot <- plot + scale_x_continuous(breaks = scales::pretty_breaks(n = xBreaks));
  }

  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.LINE_PLOT_NAME, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
  
  return(plot);
}


#########################################################
#' Creates histograms for the ranks.
#'
#' @param rankCounts1 {vector} the values for the first plot
#' @param rankCounts2 {vector} the values for the second plot
#' @param rankCounts3 {vector} the values for the third plot
#' @param rankCounts4 {vector} the values for the fourth plot
#' @param titles {list} the titles for the different plots
#' @param colors {vector} vector of colors
#'
#' @export
Plotter.createRanksHistograms <- function(rankCounts1, rankCounts2, rankCounts3, rankCounts4, titles, colors) {
  finalColors <- rep(colors, length.out = 4);
  
  plot1 <- Plotter.__createRanksHistogram(rankCounts1, titles[[1]], finalColors[[1]]);
  plot2 <- Plotter.__createRanksHistogram(rankCounts2, titles[[2]], finalColors[[2]]);
  plot3 <- Plotter.__createRanksHistogram(rankCounts3, titles[[3]], finalColors[[3]]);
  plot4 <- Plotter.__createRanksHistogram(rankCounts4, titles[[4]], finalColors[[4]]);
  
  combinedPlot <- gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2);
  print(combinedPlot);
}


#########################################################
#' Creates a histogram for the ranks.
#'
#' @param rankCounts {vector} the values to plot from the lowest to the highest
#' @param title {string} the title for the histogram
#' @param color {string} the string for the fill of a bar
#' 
#' @return {ggplot} the plot
Plotter.__createRanksHistogram <- function(rankCounts, title, color) {
  xData <- 1:length(rankCounts);
  yData <- rankCounts;
  
  return(Plotter.plotGeneralBarPlot(xData, yData, title, 
                                    Strings.RANKS, Strings.FREQUENCY, 
                                    FALSE, FALSE, TRUE, c(0, 180), TRUE, color));
}


#########################################################
#' Plots a custom line plot.
#'
#' @param xData {vector} the vector with x-axis data
#' @param yData {vector} the vector with y-axis data
#' @param group {vector} the vector with group data
#' @param plotTitle {string} the title of the plot
#' @param xAxisLabel {string} the label for the x-axis
#' @param yAxisLabel {string} the label for the y-axis
#' @param groupLabel {string} the label for the group
#' @param print {logical} show the plot or not
#' @param save {logical} save the plot as or not
#' @param points {logical} tells if points should be plotted or not
#' @param noLegend {logical} tells if no legend should be plotted
#' @param secondAxis {logical} tells if second x-Axis should be created
#'
#' @export
Plotter.plotMultipleLinePlots <- function(xData, yData, group, plotTitle, 
                                          xAxisLabel, yAxisLabel, groupLabel, 
                                          print, save, labels, colors, 
                                          points = FALSE, yLimits = vector(),
                                          noLegend = TRUE, secondAxis = FALSE) {
  # encode data
  data <- data.frame(x = xData, y = yData, g = group);
  
  # create plot
  plot <- ggplot(data, aes(x = x, y = y, color = g)) + 
    geom_line() +
    labs(
      title = plotTitle,
      x = xAxisLabel,
      y = yAxisLabel
    ) + 
    scale_color_manual(name = groupLabel, values = colors);

  # add points
  if (points) {
    plot <- plot + geom_point();
  }
  
  if (length(yLimits) > 0) {
    plot <- plot + coord_cartesian(ylim = yLimits);
  }
  
  if (noLegend) {
    plot <- plot + theme(legend.position = "none");
  }
  
  if (secondAxis) {
    plot <- plot + scale_y_continuous("", sec.axis = sec_axis(~.+10, name = derive()))
  }
  
  if (print) print(plot);
  if (save) Storer.savePlot(plot, paste(Defaults.LINE_PLOT_NAME, 
                                        Extensions.PNG, sep = Symbols.EMPTY));
}


#########################################################
#' Plots a log-normal distribution.
#'
#' @param scores {vector} the scores used for the histogram
#' @param fit {list} the to the scores fitted parameters of the distribution
#' @param color {string} the color for the density curve
#' @param print {logical} tells if the plot has to be printed or not
#'
#' @return {list} the ggplot
Plotter.getLogNormalDistributionPlot <- function(scores, fit, color, print) {
  data <- data.frame(scores = scores)
  plot <- ggplot(data, aes(x = scores)) + 
    geom_histogram(aes(y = ..density..), 
                   colour = Colors.WHITE, 
                   fill = Colors.GREEN,
                   binwidth = 1) +
    stat_function(fun = dlnorm, args = list(meanlog = fit$estimate[1], sdlog = fit$estimate[2]),
                  color = color) + 
    labs(
      x = Strings.SCORES,
      y = Strings.PROBABILITY
  );
  
  if (print) print(plot);
  
  return(plot);
}