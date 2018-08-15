#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Starts the whole program.

# libraries
library(ggplot2);  # to draw nice histograms and graphs
library(gridExtra);  # to create graphs side by side
library(fitdistrplus);  # allows fitting distributions
library(plyr);  # to split, combine and apply on data
library(reshape2);  # to transform data easily for ggplot
library(rlist);  # to get different list functions
library(scales);  # to get special axes in ggplot
library(stringr);   # to modify strings
library(VennDiagram);  # to draw Venn diagrams

source("libraries/mica-functions.R");

# implementations
if(!exists("ANALYZER_IMPORTED")) source("maths/Analyzer.R");
if(!exists("CLUSTERING_IMPORTED")) source("maths/Clustering.R");
if(!exists("CURVE_MINER_IMPORTED")) source("maths/CurvesMiner.R");
if(!exists("DATA_ANALYSIS_IMPORTED")) source("experiments/DataAnalysis.R");
if(!exists("DEFAULTS_IMPORTED")) source("system/Defaults.R");
if(!exists("ENCODER_IMPORTED")) source("visuals/Encoder.R");
if(!exists("EXPERIMENT_0_IMPORTED")) source("experiments/Experiment0.R");
if(!exists("EXPERIMENT_1_IMPORTED")) source("experiments/Experiment1.R");
if(!exists("EXPERIMENT_2_IMPORTED")) source("experiments/Experiment2.R");
if(!exists("EXPERIMENT_3_IMPORTED")) source("experiments/Experiment3.R");
if(!exists("EXPERIMENT_4_IMPORTED")) source("experiments/Experiment4.R");
if(!exists("EXPERIMENT_5_IMPORTED")) source("experiments/Experiment5.R");
if(!exists("INTERFACE_IMPORTED")) source("system/Interface.R");
if(!exists("INTERPRETER_IMPORTED")) source("system/Interpreter.R");
if(!exists("LOADER_IMPORTED")) source("system/Loader.R");
if(!exists("MATH_IMPORTED")) source("maths/Math.R");
if(!exists("PLOTTER_IMPORTED")) source("visuals/Plotter.R");
if(!exists("PRESENTATION_IMPORTED")) source("experiments/Presentation.R");
if(!exists("STORER_IMPORTED")) source("system/Storer.R");

#########################################################
#' Starts the program.
#' hint: under "Files" (File Explorer) set the files-folder 
#' with this "Main"-class contained as working directory
Main.main <- function() {
	Main.initialize();
	
  # scripts with the tasks used to create the corresponding thesis
  # Main.experiments();
  
  # conversion into the right interface format
  # Main.conversion();
  
  # interface for approaches (for use in dendrochronology chair)
  Main.interface();
  
  # for the colloquium 
  # Presentation.start();
  
  print(Strings.DONE);
}


#########################################################
#' Starts the experiment-scripts which contain the tasks which had to be executed during the thesis.
Main.experiments <- function() {
  # tests (with artificial set)
  # Experiment0.start(Paths.ARTIFICAL_1_LOW_RES_DATASET,
  #                   Defaults.CURVE_LOW_RES_NAME);
  
  # data analysis
  # DataAnalysis.start();
  
  # consensus approach
  # Experiment1.start();
  
  # bucket approach
  # Experiment2.start();
  
  # voting based approach
  # Experiment3.start();
  
  # two-step approach
  # Experiment4.start();
  
  # per tree approach
  # Experiment5.start();
}


#########################################################
#' Reconverts data for the interfaces, since the interfaces do not use Dr. Martin's and Alex computed files
#' because ring-widths and consensi were stored in seperate files (working examples - that's why it is not read from Defaults!).
Main.conversion <- function() {
  # Storer.storeInDefaultFormat(pathProfiles = "input/conversions/pass_1/",
  #                             pathCharacteristics = "input/conversions/pass_1/",
  #                             profileName = "profiles",
  #                             widthName = "ring_widths",
  #                             outputFileName = "consensus");

  # Storer.storeAllInDefaultFormat(pathProfiles = "input/conversions/pass_1/samples_profiles/",
  #                                pathCharacteristics = "input/conversions/pass_1/samples_ring_widths/");
  
  # Storer.storeAllInDefaultFormat(pathProfiles = "input/conversions/consensi_profiles/",
  #                                pathCharacteristics = "input/conversions/consensi_ring_widths/",
  #                                isSample = FALSE);  # isSample set to false to avoid a lost of the years
  
  # Storer.removeSampleYearsFromAll(consensiPath = "input/conversions/consensi/",
  #                                 samplesPath = "input/interface/pass_1/samples/");
} 


#########################################################
#' Interface for the  different approaches 
#' (working examples - that's why it is not read from Defaults!).
Main.interface <- function() {
  # Hint: Same samples dated with different approaches.
  
  # dates <- Interface.computeDatesConsensusApproach(consensusPath = "input/interface/pass_1/",
  #                                                  consensusName = "consensus",
  #                                                  samplesPath = "input/interface/pass_1/samples/",
  #                                                  scoreType = "d",
  #                                                  bestYearsMax = 5,
  #                                                  save = TRUE,
  #                                                  fileName = "datesConsensus");
  
  # dates <- Interface.computeDatesPerTreeApproach(perTreePath = "input/interface/pass_1/buckets_consensi/",
  #                                                samplesPath = "input/interface/pass_1/samples/",
  #                                                scoreType = "d",
  #                                                bestYearsMax = 5,
  #                                                save = TRUE,
  #                                                fileName = "datesPerTree");
  # 
  # dates <- Interface.computeDatesBucketApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #                                               samplesPath = "input/interface/pass_1/samples/",
  #                                               scoreType = "d",
  #                                               innerFunc = min,
  #                                               outerFunc = sum,
  #                                               bestYearsMax = 5,
  #                                               qualityMeasure = "ps",
  #                                               save = TRUE,
  #                                               fileName = "datesBucket");

  # dates <- Interface.computeDatesVotingApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #                                               samplesPath = "input/interface/pass_1/samples/",
  #                                               scoreType = "d",
  #                                               topYearsCount = 1,
  #                                               approach = "p",  # activates powerset approach
  #                                               minimumLength = 8,
  #                                               bestYearsMax = 5,
  #                                               save = TRUE,
  #                                               fileName = "datesVoting");
  
  # dates <- Interface.computeDatesTwoStepApproach(bucketsPath = "input/interface/pass_1/buckets_consensi/",
  #                                                samplesPath = "input/interface/pass_1/samples/",
  #                                                scoreTypeCharacteristic = "p",
  #                                                scoreTypeBucket = "d",
  #                                                topYearsCount = 20,
  #                                                bestYearsMax = 5,
  #                                                qualityMeasures = "s",
  #                                                save = TRUE,
  #                                                fileName = "datesTwoStep");
}


#########################################################
#' Does all the initializations.
Main.initialize <- function() {
  initMica("libraries");  # last line in mica script for initialization has to be removed
}


# start program
Main.main();
