#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

#' Executes Unit-Tests.

# libraries
library(ggplot2);  # to draw nice histograms and graphs
library(gridExtra);  # to create graphs side by side
library(plyr);  # to split, combine and apply on data
library(reshape2);  # to transform data easily for ggplot
library(rlist);  # to get different list functions
library(scales);  # to get percentages in plots
library(stringr);   # to modify strings
library(testthat);  # to do the unit-test
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
if(!exists("LOADER_IMPORTED")) source("system/Loader.R");
if(!exists("MATH_IMPORTED")) source("maths/Math.R");
if(!exists("PLOTTER_IMPORTED")) source("visuals/Plotter.R");
if(!exists("PRESENTATION_IMPORTED")) source("experiments/Presentation.R");
if(!exists("STORER_IMPORTED")) source("system/Storer.R");

#' hint: under "Files" (File Explorer) set the files-folder above as working directory (where the "Main"-class contained)
initMica("libraries");
test_results <- test_dir("tests", reporter="summary")