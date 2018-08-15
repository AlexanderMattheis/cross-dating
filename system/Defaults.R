#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

DEFAULTS_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Used to store all constants and default parameters.
#' Functional strings like "dotted" or "histogram" aren't stored here,
#' because it would decrease the speed in the whole development.
#' Hint: Alphabetically and Numerically sorted!

# approaches
Approaches <- c("Ring-Width", "Consensus", "Two-Step", "Voting", "Per-Tree", "Bucket");

# check
Check.NUMBER_OF_CURVES <- "#Curves:";
Check.NUMBER_OF_SAMPLES <- "#Samples:";
Check.SAMPLE_LENGTH <- "Sample length:";

# colors
Colors.BLACK <- "#000000";
Colors.BLUE <- "#0bbfd8";
Colors.BROWN <- "#ba6526";
Colors.GRAY <- "azure3";
Colors.GREEN <- "#a0bf1c";
Colors.ORANGE <- "#ffb300";
Colors.PINK <- "#ff007f";
Colors.RED <- "#e63900";
Colors.VIOLET <- "#be0086";
Colors.WHITE <- "white";
Colors.YELLOW <- "#ffe22f";

  # dark colors
Colors.DARK_BLUE <- "#076195";
Colors.DARK_BROWN <- "#7d6526";
Colors.DARK_GREEN <- "#91af18";
Colors.DARK_RED <- "#ac2e00";
Colors.DARK_VIOLET <- "#862d86";

  # light colors
Colors.LIGHT_BLUE <- "#0b8eda";

Palettes.BLUE <- c("#0b8eda", "#0ba9da", Colors.BLUE, Colors.LIGHT_BLUE);
Palettes.OTHER_NATURE <- c("#e70073", "#ffb300", "#ff6666", Colors.YELLOW);
Palettes.NATURE <- c(Colors.GREEN, Colors.RED, Colors.BLUE, Colors.VIOLET,
                     Colors.ORANGE, Colors.DARK_GREEN, "#cd3700", Colors.DARK_BLUE,
                     Colors.DARK_VIOLET); 
Palettes.RED <- c("#e60000", "#e70073", "#e4b0cA");
Palettes.VIOLET <- c("#862d86", "#982074", Colors.VIOLET, "#d40094"); 

Palettes.NONE <- list();

# defaults (MICA encoding etc.) 
# HINT: It is dangerous to do changes on the defaults! The whole program can be broken!
Defaults.ALL_VOTED_FOR_THE_SAME <- -1;
Defaults.BAR_PLOT_NAME <- "bar_plot";
Defaults.BACK_TO_BACK_HISTOGRAM <- "Back-To-Back";
Defaults.BOX_PLOT_NAME <- "box_plot";
Defaults.CONSENSUS_CURVE_OSTALB_FILENAME <- "PCAB_Ostalb_GD-pruned_Consensus.csv";
Defaults.CONSENSUS_PLOT_PDF_NAME <- "consensi";
Defaults.CORRELATION_GLEICHLAUF <- "gleichlauf";
Defaults.CORRELATION_KENDALL <- "kendall";
Defaults.CORRELATION_PEARSON <- "pearson";
Defaults.CORRELATION_SPEARMAN <- "spearman";
Defaults.CORRELATION_T_VALUE <- "t-value";
Defaults.CURVE_HIGH_RES_NAME <- "curveHighRes";  # generic name
Defaults.CURVE_LOW_RES_NAME <- "curveLowRes";  # generic name
Defaults.CURVE_MID_RES_NAME <- "curveMidRes";  # generic name
Defaults.CURVE_OSTALB_END_YEAR <- 2004;  # default: 2004 in Ostalb dataset
Defaults.CURVE_OSTALB_NAME <- "_MICA-cons";  # generic name
Defaults.CURVE_OSTALB_START_YEAR <- 1916;  # default: 1916 in Ostalb dataset
Defaults.CURVE_MEAN_ABSOLUTE_DISTANCE <- 2;
Defaults.CURVE_RMSD <- 0;
Defaults.CURVE_SLOPE_RMSD <- 1;
Defaults.CURVE_SLOPE_MEAN_ABSOLUTE_DISTANCE <- 3;
Defaults.DELTA_SCORE <- expression(Delta~"Score");
Defaults.DENSITY_NAME <- "density";
Defaults.DISTRIBUTION_GAMMA <- "gamma";
Defaults.DISTRIBUTION_LOGNORMAL <- "lognormal";
Defaults.DISTRIBUTION_NORMAL <- "normal";
Defaults.DISTRIBUTION_WEIBULL <- "weibull";
Defaults.DOUBLE_WEIGHTING_COMMAND <- "d";
Defaults.GENERIC_DENSITY_COLUMN_NAME <- "GD";
Defaults.GENERIC_WIDTH_COLUMN_NAME <- "width";
Defaults.HISTOGRAM <- "histogram";
Defaults.INTERPOLATION_SIZE <- 250;
Defaults.LINE_PLOT_NAME <- "line_plot";
Defaults.LOG_WEIGHTING_COMMAND <- "l";
Defaults.MASTER_CHRONOLOGY_OSTALB_FILENAME <- "PCAB_Ostalb_GD-pruned_Consensus-sample-free.csv";
Defaults.NUMBER_OF_DISTANCE_FUNCTIONS <- 4;
Defaults.OUTLIER <- -10;
Defaults.PART_NAME <- "part";
Defaults.POWER_SET_APPROACH_COMMAND <- "p";
Defaults.P_VALUE_GROUP <- "length-10\nsample";
Defaults.P_VALUE_NAME <- "pValue";
Defaults.P_VALUES <- "p-Value";
Defaults.P_VALUES_COMMAND <- "p";
Defaults.RANK_NAME <- "rank";
Defaults.ROC_PLOT_NAME <- "roc_plot";
Defaults.ROUNDING_DIGITS <- 4;  # tells on how many digits scores should be rounded
Defaults.ROUNDING_DIGITS_2 <- 8;  # tells on how many digits p-values should be rounded
Defaults.SAMPLE_NAME <- "sample";
Defaults.SCALING <- 1;  # set back to 1, because there were no real influence noticeable
Defaults.SCATTER_PLOT_NAME <- "scatter_plot";
Defaults.SCORES_COMMAND <- "s";
Defaults.SCORE_NAME <- "score";
Defaults.TASK_A <- "a";
Defaults.TASK_B <- "b";
Defaults.TASK_C <- "c";
Defaults.TASK_D <- "d";
Defaults.TASK_E <- "e";
Defaults.TASK_F <- "f";
Defaults.TASK_G <- "g";
Defaults.TASK_H <- "h";
Defaults.TASK_KENDALL_CHAR <- "t";
Defaults.TASK_GLEICHLAUF_CHAR <- "g";
Defaults.TASK_PEARSON_CHAR <- "p";
Defaults.TASK_RING_WIDTH_BUCKET <- "w";
Defaults.TASK_SPEARMAN_CHAR <- "r";
Defaults.TASK_T_VALUE_CHAR <- "v";
Defaults.TUKEYS_C <- 9;
Defaults.VIOLINE_PLOT_NAME <- "violin_plot";
Defaults.YEAR_NAME <- "year";

# exceptions
Exceptions.CURVE_YEARS_NOT_IDENTICAL <- "The curve years of both curves have to be identical!";
Exceptions.CURVES_NUM_PROFILES <- "The files storing the num of profiles have different lengths!";
Exceptions.TOP_YEARS <- "Minimum number of top-years is two!";
Exceptions.WRONG_EXTENSION <- "One or more of your files have a wrong extension!";

# experiments
Experiment0.CONSENSUS_OF_TREE <- list(1,1,1,1,
                                      2,2,2,2,
                                      3,3,3,3,
                                      4,4,4,4,
                                      5,5,5,5);

Experiment0.DISTANCE_FUNCTIONS <- list(0, 1, 2, 3);

Experiment0.PATTERNS <- list(c(5, 16), c(16, 25), c(25, 35), c(33, 46),
                             c(4, 12), c(17, 25), c(26, 40), c(12, 18),
                             c(25, 39), c(26, 36), c(58, 70), c(32, 42),
                             c(19, 29), c(17, 27), c(63, 72), c(39, 45),
                             c(52, 60), c(65, 72), c(25, 37), c(10, 20));

Experiment0.PATTERN_OF_TREE <- list("0111", "0121", "0131", "0141", 
                                    "0211", "0221", "0231", "0281",
                                    "0321", "0351", "0361", "0381",
                                    "0411", "0421", "0431", "0451",
                                    "0521", "0541", "0551", "0561");

# expressions
Expression.START <- "^";
Expression.NUMBERS <- "[0-9]+"; 
Expression.YEARS <- "[0-9]{4}-[0-9]{4}$";
Expression.YEARS_IN_BETWEEN <- "[0-9]{4}-[0-9]{4}";

# extensions
Extensions.CSV <- ".csv";
Extensions.META <- ".meta";
Extensions.PDF <- ".pdf";
Extensions.PNG <- ".png";
Extensions.TIFF <- ".tiff";

# files
Files.INTERSECTION <- "Venn_Max_Intersection_";
Files.OSTALB_CHRONOLOGY <- "PCAB_Ostalb_GD-pruned_Consensus-sample-free";
Files.OSTALB_CONSENSUS <- "PCAB_Ostalb_GD-pruned_Consensus";  # the same with samples
Files.OSTALB_CONSENSUS_NUM_PROFILES <- "PCAB_Ostalb_GD-NumProfiles-Consensus"; 
Files.OSTALB_CONSENSUS_SAMPLE_FREE_NUM_PROFILES <- "PCAB_Ostalb_GD-NumProfiles-Consensus-sample-free";
Files.RANK_HISTOGRAM <- "Rank_Histogram";
Files.SCORES_FOLDER_NAME <- "PCAB_Ostalb_GD-pruned_Scores_";
Files.SINGLE_CONSENSUS_PROFILE <- "profile_Consensus-sample-free";
Files.SINGLE_SAMPLE_PROFILE <- "profile_sample-1939-1948";
Files.VENN_DIAGRAM <- "Venn_Diagramm_";

# math
Maths.INFINITY <- 2^30;
Maths.EPSILON <- 0.0001; 

# paths
Paths.ALL_VS_ALL_DISTANCES <- "input/all_vs_all_distances/";
Paths.ARTIFICAL_1_HIGH_RES_DATASET <- "input/artificial/artifical1_high_res/";
Paths.ARTIFICAL_1_LOW_RES_DATASET <- "input/artificial/artifical1_low_res/";
Paths.ARTIFICAL_1_MID_RES_DATASET <- "input/artificial/artifical1_mid_res/";
Paths.BUCKETS <- "per_year_profiles/";
Paths.INPUT <- "input/";
Paths.PASSES <- "input/passes/pass_";  # here the consensi are stored
Paths.PASSES_BUCKET_MIN <- "input/min_bucket_passes/pass_";
Paths.PASSES_BUCKET_MIN_LENGTH_1 <- "input/min_bucket_passes_lenght_1/pass_";
Paths.PASSES_BUCKET_MIN_LENGTH_15 <-  "input/min_bucket_passes_length_15/pass_";
Paths.PASSES_BUCKET_MIN_LENGTH_5 <- "input/min_bucket_passes_length_5/pass_";  # here the consensi are stored
Paths.PASSES_CLUSTERED <-  "input/clustered/pass_";
Paths.PASSES_COUNTS <- "input/passes_bucket_counts/pass_";
Paths.PASSES_FILTERED <- "input/filtered_bucket_passes/pass_";
Paths.PASSES_MIN_MIN <- "input/min_min_passes/pass_";
Paths.PASSES_MAX_MIN <- "input/max_min_passes/pass_";
Paths.PASSES_MAX_DENSITY <- "input/passes_max_density/pass_";
Paths.PASSES_MAX_DENSITY_LENGTH_15 <- "input/passes_max_density_length_15/pass_";
Paths.PASSES_MAX_DENSITY_LENGTH_5 <- "input/passes_max_density_length_5/pass_";
Paths.PASSES_LENGTH_1 <- "input/passes_length_1/pass_";
Paths.PASSES_LENGTH_15 <- "input/passes_length_15/pass_";
Paths.PASSES_LENGTH_5 <- "input/passes_length_5/pass_";
Paths.PASSES_PER_TREE <- "input/passes_per_tree/pass_"; 
Paths.PASSES_PER_TREE_MAX_DENSITY <- "input/passes_per_tree_max_density/pass_";
Paths.PASSES_PER_TREE_MAX_DENSITY_LENGTH_15 <- "input/passes_per_tree_max_density_length_15/pass_";
Paths.PASSES_PER_TREE_LENGTH_15 <- "input/passes_per_tree_length_15/pass_";
Paths.PASSES_PER_TREE_LENGTH_5 <- "input/passes_per_tree_length_5/pass_";
Paths.PASSES_PER_TREE_SUBDIV_15_LENGTH_15 <- "input/passes_per_tree_subdivisions_15_length_15/pass_";
Paths.PASSES_PER_TREE_SUBDIV_2 <- "input/passes_per_tree_subdivisions_2/pass_";
Paths.PASSES_PER_TREE_SUBDIV_3_LENGTH_15 <- "input/passes_per_tree_subdivisions_3_length_15/pass_";
Paths.PASSES_PER_TREE_SUBDIV_5_LENGTH_15 <- "input/passes_per_tree_subdivisions_5_length_15/pass_";
Paths.PASSES_PER_TREE_SUBDIV_5 <- "input/passes_per_tree_subdivisions_5/pass_";
Paths.PASSES_RING_WIDTH <- "input/passes_ring_width/pass_"; 
Paths.PASSES_RING_WIDTH_BUCKET <- "input/passes_ring_width_bucket/pass_"; 
Paths.PASSES_RING_WIDTH_BUCKET_15 <- "input/passes_ring_width_bucket_length_15/pass_"; 
Paths.PASSES_RING_WIDTH_BUCKET_5 <- "input/passes_ring_width_bucket_length_5/pass_"; 
Paths.PASSES_RING_WIDTH_LENGTH_15 <- "input/passes_ring_width_length_15/pass_";
Paths.PASSES_RING_WIDTH_LENGTH_5 <- "input/passes_ring_width_length_5/pass_";
Paths.PASSES_RING_WIDTH_PER_TREE <- "input/passes_ring_width_per_tree/pass_";
Paths.PASSES_VOTING <- "input/voting_passes/pass_";
Paths.PASSES_VOTING_LENGTH_15 <- "input/voting_passes_length_15/pass_";
Paths.PASSES_VOTING_LENGTH_5 <- "input/voting_passes_length_5/pass_";
Paths.OSTALB_DATASET <- "input/PCAB_Ostalb_GD-pruned_MICA-per-tree/";
Paths.OSTALB_DATASET_WIDTHS <- "input/PCAB_Ostalb_GD-ring_widths/";
Paths.OSTALB_DATASET_MAXIMUM_DENSITIES <- "input/PCAB_Ostalb_GD-max_densities/";
Paths.OSTALB_TESTSET_FOLDER <- "PCAB_Ostalb_GD-pruned_Test-samples/";
Paths.OUTPUT <- "output/";
Paths.PATH <- "/";
Paths.TEST_CURVES <- "../input/tests/test_curves/";
Paths.TEST_DATASET <- "../input/tests/PCAB_Ostalb_GD-pruned_MICA-per-tree/";
Paths.TEST_PASSES_PER_TREE <- "../input/tests/passes_per_tree/pass_"; 
Paths.TEST_SCORES <- "../input/tests/test_scores/";
Paths.TEST_SCORES_2 <- "../input/tests/test_scores_2/";

# parameters (which can be set differently)
Parameters.LENGTH_TO_REMOVE <- 9;
Parameters.LIMIT_PROFILES_PER_CONSENSUS <- 25;
Parameters.RANKS_TO_OUTPUT <- 5;  # for how many ranks the number of samples should be displayed in the violine plots
Parameters.SCORE_TYPES_TO_LOOK_AT <- c(a = "a", b = "b", c = "c", d = "d");
Parameters.SELECTION_NUMBER <- 12;  # tells how many bad and good curves should be selected
Parameters.TEST_SET_SIZE <- 20;  # dependant on Experiment 0 parameters
Parameters.TOP_RANK <- 5;
Parameters.TRIES_TO_FIND_A_SAMPLE <- 20;

  # values used to compute Ostalb per tree dataset as comment
  # descriptions taken from https://github.com/BackofenLab/MICA
Parameters.DIST_FUNCTION <- 3;  # 3
Parameters.DIST_SAMPLE <- 1000;  # 500, number of samples used to compute distance (and so the consensus)
Parameters.DIST_WARP_SCALING <- 0;  # 0, if > 0, then multiplied with distance
Parameters.MAX_WARPING_FACTOR <- 3;  # 2.5, maximally allowed interval distortion
Parameters.MAX_REL_X_SHIFT <- 0.1;  # 0.1, maximally allowed relative shift of x-coordinates
Parameters.MIN_PROFILES_PER_YEAR <- 3;  # 3, minimal number or profiles used for consensus of a tree
Parameters.MIN_REL_INTERVAL_LENGTH <- 0.06;  # 0.05, minimal relative interval length which is still considered for decomposition 
Parameters.MIN_REL_MIN_MAX_DIST <- 0.04;  # 0.02, minimal distance between neighbored minima and maxima
Parameters.MIN_REL_SLOPE_HEIGHT <- 0.04;  # 0.02, minimal relative slope value for inflection points
Parameters.REFERENCE <- 0;  # 0, index of a curve
Parameters.OUT_SLOPE <- FALSE;  # FALSE, output or not output computed slope values

Parameters.RANKING_DIST_FUNCTION <- 3;
Parameters.RANKING_DIST_SAMPLE <- 250;
Parameters.RANKING_DIST_WARP_SCALING <- 0;
Parameters.RANKING_MAX_WARPING_FACTOR <- 3;
Parameters.RANKING_MAX_REL_X_SHIFT <- 0.1;
Parameters.RANKING_MIN_REL_INTERVAL_LENGTH <- 0.06;
Parameters.RANKING_MIN_REL_MIN_MAX_DIST <- 0.04;
Parameters.RANKING_MIN_REL_SLOPE_HEIGHT <- 0.04;
Parameters.RANKING_REFERENCE <- 1;
Parameters.RANKING_OUT_SLOPE <- FALSE;

# strings
Strings.APPROACH <- "Approach";
Strings.BEST_RANKED <- "Best ranked\nsamples";
Strings.BUCKET_BASED <- "bucket-based";
Strings.BUG <- "Error! Your dataset contains same profiles twice.";
Strings.CORRELATION_COEFFICIENT <- "Pearson Correlation\nCoefficient";
Strings.CONSENSI <- "Consensi";
Strings.CONSENSUS_1_NAME <- "Consensus";
Strings.CONSENSUS_2_NAME <- "Consensus\nsample-free";
Strings.COUNT <- "Count:    ";
Strings.COUNT_2 <- "Count";
Strings.COUNT_LOWER_CASE <- "count";
Strings.DATASET_CORRUPT <- "Dataset contains same profiles twice!";
Strings.DENSITY <- "Density (GD)";
Strings.DELTA_PEAK <- expression(Delta~"Peak");
Strings.DELTA_PREDICTION_ERROR <- expression(Delta~"Prediction Error");
Strings.DISTANCES <- "Distance";
Strings.DONE <- "DONE!";
Strings.FALSE_POSITIVE_RATE <- "fallout(d)";
Strings.FINISHED <- "Finished!";
Strings.FREQUENCY <- "Frequency";
Strings.IS_RANK_1 <- "Is rank 1?";
Strings.LENGTH <- "Series-Length";
Strings.MEAN <- "Mean:    ";
Strings.MEANS <- "Truth Means";
Strings.MEAN_VALUE <- "Mean";
Strings.MEDIAN <- "Median:      ";
Strings.NOT_IMPLEMENTED_YET <- "Not implemented yet!";
Strings.OUTLIERS <- "Outcasts:   "; 
Strings.NONE <- "[None]";
Strings.NUMBER <- "Number";
Strings.NUMBER_OF_TREES <- "Number of Trees";
Strings.NUMBER_OF_SUBSAMPLES <- "Number of Subsamples";
Strings.NUMBER_RANK_1 <- "Number of rank 1";
Strings.NUM_RANK_ <- "#Rank ";
Strings.NUM_RANK_BIGGER <- "#BiggerRanks ";
Strings.PEAK_MEAN <- "Peaks";
Strings.PEAK_RANK <- "Peak-Rank";
Strings.POINT <- "Point";
Strings.PREDICTION_ERROR <- "Prediction Error";
Strings.PREDICTED_CORRECT_YEAR_RANKS <- "Predicted Rank\nof Start Year";
Strings.PREDICTED_RANK <- "Rank";
Strings.PROBABILITY <- "Probability";
Strings.PROFILES <- "Sample-Profiles:";
Strings.P_VALUES <- "p-Values";
Strings.QUALITY <- "Mean Predicted\nStart Year Rank";
Strings.QUALITY_2 <- "Correct Rated Samples";
Strings.QUALITY_3 <- "Top 5 Rated Samples";
Strings.RANK <- "Rank";
Strings.RANKS <- "Ranks";
Strings.RANK_1 <- "rank 1";
Strings.RANK_2 <- "rank 2";
Strings.RATIO <- "Ratio";
Strings.RELATIVE_QUALITY <- "Relative Quality";
Strings.REPLACEMENT <- " / ";
Strings.REPLACEMENT_2 <- " - ";
Strings.REPLACEMENT_3 <- "min_bucket_passes";
Strings.REPLACEMENT_4 <- "Scores_d";
Strings.REPLACEMENT_5 <- "Score_d";
Strings.RING_WIDTH_BASED <- "ring-width based";
Strings.RUNTIME <- "Runtime (seconds)";
Strings.RUNTIMES <- "Runtimes";
Strings.SAMPLES <- "Samples";
Strings.SCORE <- "Score";
Strings.SCORES <- "Scores";
Strings.SCORES_MEAN <- "Scores";
Strings.SCORE_LOWER_CASE <- "score";
Strings.SEPERATION <- "----";
Strings.TIME <- "Relative X";
Strings.TO_REPLACE_1 <- "Score_d.csv";
Strings.TO_REPLACE_2 <- "_MICA-cons";
Strings.TO_REPLACE_3 <- "passes_ring_width_bucket";
Strings.TO_REPLACE_4 <- "Scores_p";
Strings.TO_REPLACE_5 <- "Score_p";
Strings.TO_REPLACE_6 <- "passes_max_density";
Strings.TO_REPLACE_7 <- "Scores_v";
Strings.TO_REPLACE_8 <- "Score_v";
Strings.TRUE_POSITIVE_RATE <- "recall(d)";
Strings.TRUTH <- "Truth";
Strings.TYPE <- "Type";
Strings.VARIANCE <- "Var:   ";
Strings.YEARS <- "Year";
Strings.YEAR_LOWER_CASE <- "year";
Strings.WIDTH <- "Ring-Width";
Strings.WINDOW_LENGTH <- "Window Length";
Strings.WORST_RANKED <- "Worst ranked\nsamples";
Strings.ZERO_SCORE <- "Zero Score!";

Titles.APPROACH <- "Approach";
Titles.CHARACTERISTIC <- "Non-characteristic and below characteristic statistics";
Titles.CHARACTERISTIC_3 <- "Correct Dating in Best 3 Rankings";
Titles.CHARACTERISTIC_5 <- "Correct Dating in Best 5 Rankings";
Titles.CHARACTERISTIC_10 <- "Correct Dating in Best 10 Rankings";
Titles.YEARS_IN_MAXIMUM_INTERSECTION <- "Years in Maximum\nIntersection-Subset from";
Titles.COMPUTING <- "Computing single derivate/normalization etc.";
Titles.CONSECUTIVE <- "Years consecutive in file?";
Titles.CONSENSUS <- "Consensus";
Titles.CURRENTLY_COMPUTED <- "Printing currently preprocessed sample index.";
Titles.CURRENTLY_PROCESSED <- "Currently processed table:";
Titles.CURRENT_PROCESSED_SAMPLE <- "Printing currently processed sample:";
Titles.CURRENTLY_PROCESSED_YEAR <- "Printing currently processed year:";
Titles.CURVES <- "Curves";
Titles.DELTA_PREDICTION_ERROR <- expression(Delta~"Prediction-Error for"~Delta~"Scores");
Titles.DIFFERENCES_IN_YEARS <- "Printing Years in which number(widths) != number(profiles):";
Titles.DISTRIBUTION_OF_YEARS <- "Distribution of Years in Max. Intersection-Set (10 Worst / 10 Best)";
Titles.ENHANCEMENT <- "Ranks Before and After";
Titles.EQUIDISTANT <- "Values equidistant in file?";
Titles.INCREMENTED_YEARS_RATIO <- "Ratio of Incremented Years";
Titles.LOADING_FROM_DISK <- "Loading files from disk..";
Titles.MEANS_IN_RANGES <- "Means in Ranges of"; 
Titles.MEANS_IN_RANGES_25 <- "Means in Ranges of 25";
Titles.MEANS_IN_RANGES_0_89 <- "Means in Ranges of 8/9";
Titles.MEANS_IN_RANGES_0.009 <- "Means in Ranges of 0.009";
Titles.MEANS_IN_RANGES_11_900 <- "Means in Ranges of 11/900";
Titles.METHOD <- "Method";
Titles.NEIGHBOURED_SCORES <- "Neigboured Profiles' Scores";
Titles.NUMBER_OF_POINTS <- "Printing number of points per year:";
Titles.NUMBER_OF_PROFILES <- "Printing number of profiles per year:";
Titles.NUMBER_OF_PROFILES_FINAL <- "Total Number of profiles:";
Titles.PREDICTION_ERROR <- expression("Prediction-Error for"~Delta~"Scores");
Titles.PROFILES_AND_YEARS <- "Number of profiles and years";
Titles.PROFILES_DISTANCES_PER_YEAR <- "Distances between Normalized Profiles";
Titles.PROFILES_MEAN_DISTANCES_PER_YEAR <- "Mean-Distances between Normalized Profiles";
Titles.PROFILES_PER_YEAR <- "Profiles per Year";
Titles.PROFILES_PER_YEAR_AFTER <- "After Clustering";
Titles.PROFILES_PER_YEAR_BEFORE <- "Before Clustering";
Titles.P_VALUES_BELOW_1_PERCENT <- "p-Values below 1%";
Titles.P_VALUES_BELOW_0_1_PERCENT <- "p-Values below 0.1%";
Titles.QUALITY_FACTOR <- "Quality Factor Peak Difference";
Titles.QUALITY_FACTOR_2 <- "Quality Factor Score Difference";
Titles.QUALITY_FACTOR_3 <- "Quality Factor P-Values";
Titles.RANKS <- "Combined Approach Ranks";
Titles.RANKS_FOR_DELTA_SCORES <- expression("Ranks for"~Delta~"Scores");
Titles.RANKS_FOR_P_VALUES <- "Ranks for p-Values";
Titles.RANKS_FREQUENCIES <- "Ranks Frequencies";
Titles.SCORE_FREQUENCY_DIAGRAM_START <- "Scores for Pattern with ";
Titles.SCORE_FREQUENCY_DIAGRAM_END <- " Profile(s)";
Titles.START_CONSENSUS_COMPUTATION <- "Started consensus computation:";
Titles.STEADY <- "Years are steady?";
Titles.TREES_PER_RANK_1 <- "Trees per Rank 1";
Titles.TRUTH_MEANS <- "Truth Means of Different Approaches";
Titles.VIOLINE_RANKS_PLOT <- "Score-Rankings with Different Measuring Methods";
Titles.VIOLINE_RANKS_PLOT_2 <- "Score-Rankings with Different Series-Lengths";
Titles.VIOLINE_RANKS_PLOT_3 <- "Histograms-Rankings with Different Measuring Methods";
Titles.WIDNOW_QUALITY <- "Window Length Dependency";
Titles.WITH_AND_WITHOUT_SAMPLES_START <- "Consensus With and Without Samples: ";
Titles.YEAR_COUNT <- "Frequency for";

MultiTitles.RANK_HISTOGRAMS <- c("Bucket-Based", "t = 1", "t = 2", "t = 3");
MultiTitles.RANK_HISTOGRAMS_2 <- c("t = 2", "t = 3", "t = 4", "t = 5");
MultiTitles.RANK_HISTOGRAMS_3 <- c("a", "b", "c", "d");

# symbols
Symbols.ARTIFICIAL_DATA_SEPARATOR <- " ";
Symbols.BRACKET_LEFT <- "(";
Symbols.BRACKET_RIGHT <- ")";
Symbols.DOT <- ".";
Symbols.EMPTY <- "";
Symbols.HYPHEN <- "-";
Symbols.NEW_LINE <- "\n";
Symbols.REAL_DATA_SEPARATOR <- ",";
Symbols.SLASH <- "/";
Symbols.SPACE <- " ";
Symbols.UNDERSCORE <- "_";