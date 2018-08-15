#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

MATH_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Contains mathematical functions for common operations like interpolation, 
#' derivates, distance and normalization functions.

#########################################################
#' Computes the discrete derivates or normalization for each sample-profile in the curvesData.
#'
#' @param curvesData {list} data of multiple curve-files (profiles with its years)
#' @param normalize {logical} tells if y-values should be normalized
#' (subtracting the mean-y and dividing by the deviation of y)
#' @param derivate {logical} this parameter can be set to compute score on slopes
#' (computed after the normalization)
#' 
#' @return {list} the list of slope-profiles with their years
#' @export
Math.computeOperation <- function(curvesData, normalize, derivate) {
  transformedProfiles <- list();
  
  # if (length(curvesData) > 1) {
  #   print(Titles.CURRENTLY_COMPUTED);
  # } else {
  #   print(Titles.COMPUTING);
  # }
  
  # iterate over all curves
  for (i in 1:length(curvesData)) {
    # print(i);
    
    # get information
    sample <- curvesData[[i]]$profiles;
    
    transformed <- list();
    
    # iterate over all profiles
    for (j in 1:length(sample)) {
      profile <- sample[[j]];
      
      if (normalize) profile <- Math.computePureValues(profile);
      if (derivate) profile <- Math.computeDerivate(profile);
      
      transformed <- rlist::list.append(transformed, profile);
    }
    
    # add information
    data <- list(profiles = transformed);
    transformedProfiles <- rlist::list.append(transformedProfiles, data);
  }
  
  return(transformedProfiles);
}


#########################################################
#' Returns a profile in which every y-value is normalized
#' by subtracting the y-mean and dividing through
#' the deviation of y.
#'
#' @param profile {list} the profile you want normalize
#'
#' @return {list} z-scores
#' @export
Math.computePureValues <- function(profile) {
  profileVector <- unlist(profile);
  
  yMean <- mean(profileVector);
  yDeviation <- sd(profileVector);  # dominator is "n-1"
  
  return(as.list((profileVector - yMean) / yDeviation));  # from each value subtract mean and divide by standard deviation
}


#########################################################
#' Computes the discrete derivate of the given profile
#' by using the lm-function in R i.e. a line
#' is computed between points (y_{i-2}, ..., y_i, ..., y_{i+2})
#' to simulate a tangent. The slope of this line is then used 
#' as gradient value s_i. 
#' Hint: Idea taken from source below.
#'
#' @param profile {list} the y-values of a discrete function
#'
#' @return {list} the values of the slope
#' @export
#' @source (idea)
#' Bender, Bela J., et al. 
#' "Microstructure alignment of wood density profiles: 
#' an approach to equalize radial differences in growth rate." 
#' Trees 26.4 (2012): 1267-127
Math.computeDerivate <- function(profile) {
  gradient <- list();
  
  # MICA based
  slopeCurve <- getSlope(getEquiX(unlist(profile)), unlist(profile));  # getEquiX because MICA xWarped was between 0 and 1
  
  # R based
  # # iterate over all y-values in the profile
  # for (i in 3:(length(profile)-2)) {
  #   x <- -2:2;
  #   y <- unlist(profile[(i-2):(i+2)]);
  # 
  #   func <- lm(y ~ x);  # approximate a line between points in profile[(i-2):(i+2)]
  # 
  #   # hint: first coefficient is the intercept with y-axis
  #   slope <- func$coefficients[[2]];  # i.e. y = intercept + slope*x = b + m*x
  #   gradient <- rlist::list.append(gradient, slope);
  # }
  
  return(slopeCurve$y);
}


#########################################################
#' Interpolates points to a profile until it has reached the
#' desired number of points.
#' Hint: The desired number of points has to be bigger
#' than the number of points in the profile.
#'
#' @param profile {list} the y-values of a profile
#' @param desiredNumberOfPoints {numerical} the desired number
#' of points
#'
#' @return {list} profile with the desired number of points
#' @export
Math.interpolatePoints <- function(profile, desiredNumberOfPoints) {
  x <- (0:(length(profile)-1))/(length(profile)-1);
  y <- unlist(profile);
  
  interpolation <- approx(x, y, n = desiredNumberOfPoints);
  
  return(as.list(interpolation$y));
}


#########################################################
#' Returns the average distance between the two profiles.
#'
#' @param firstProfile {list} the y-values of the first profile
#' @param secondProfile {list} the y-values of the second profile
#'
#' @return {numerical} the mean distance
Math.getMeanDistance <- function(firstProfile, secondProfile) {
  firstProfileVector <- unlist(firstProfile);
  secondProfileVector <- unlist(secondProfile);
  
  numberOfPoints <- length(firstProfile);
  sum <- sum(abs(firstProfileVector - secondProfileVector));
  
  return(sum/numberOfPoints);  # division is not really needed
}


#########################################################
#' Returns the Euclidean Distance.
#'
#' @param firstProfile {list} the values of the first profile
#' @param secondProfile {list} the values of the second profile
#'
#' @return {numerical} the distance between both profiles
#' @export
Math.getEuclideanDistance <- function(firstProfile, secondProfile) {
  firstProfileVector <- unlist(firstProfile);
  secondProfileVector <- unlist(secondProfile);
  
  return(dist(rbind(firstProfileVector, secondProfileVector)));
}


#########################################################
#' Returns a profile in which every y-value is normalized
#' by subtracting the y-mean and dividing through
#' the deviation of y.
#'
#' @param profiles {list} the list of profiles you want normalize
#'
#' @return {list} normalized profiles
Math.computePureProfiles <- function(profiles) {
  normalizedProfiles <- list();
  
  for (i in 1:length(profiles)) {
    normalizedProfile <- Math.computePureValues(profiles[[i]]);
    normalizedProfiles <- rlist::list.append(normalizedProfiles, normalizedProfile);
  }
  
  return(normalizedProfiles);
}


#########################################################
#' Implements Tukey's one-step biweight algorithm for robust mean calculation
#' using a w-estimator instead of a m-estimator.
#' Hint: Parameter c is by default set on value 9. 
#'
#' @param x {vector} the values 
#' for which the value has to be computed
#' @param c {numerical} tuning parameter
#'
#' @return {numerical} the function value
#' @source Statistical Algorithms Description Document, 2002, Affymetrix
#' at page 22
Math.tukeysBiweightRobustMean <- function(x, c = 9) {
  median <- median(x);
  mad <- mad(x, constant = 1);  # median(c * |x_i - median(x)|), hint: default method multiplies with constant c = 1.4826
  
  # epsilon to avoid division by zero
  uValues <- (x-median)/(c * mad + Maths.EPSILON);  # analogous to (x-mean)/sd
  
  weights <- sapply(uValues, Math.weight);  # where outlier weights removed because their values are 0

  mean <- sum(weights*x)/sum(weights);  # where outlier weights removed because their values are 0
  
  return(mean);  # compute mean
} 


#########################################################
#' Implements Tukey's weighting function.
#'
#' @param u {numerical} the value for which a weight has to be computed
#'
#' @return {numerical} a weight
Math.weight <- function(u) {
  if (abs(u) <= 1) {
    return((1-u^2)^2);
  } 
  
  return(0);
}


#########################################################
#' Computes the corrected Gleichlauf-Correlation Coefficient 
#' as proposed by Buras and Wilmking 2014/2015.
#'
#' @param x {vector} a numeric vector
#' @param y {vector} a numeric vector of same size as x
#'
#' @return {numeric} the Gleichlauf-Correlation Coefficient (value between -1 and 1)
#' @source https://doi.org/10.1016/j.dendro.2015.03.003
Math.gleichlauf <- function(x, y) {
  xDiff <- Math.__differencesBetweenNeighbouredValues(x);
  yDiff <- Math.__differencesBetweenNeighbouredValues(y);
  
  xSignReward <- sapply(xDiff, Math.__translateToReward);
  ySignReward <- sapply(yDiff, Math.__translateToReward);
  
  gleichlauf <- 1 - (sum(abs(xSignReward - ySignReward)) / length(xSignReward));  # 1 - Gegenläufigkeit
  return(gleichlauf);
}


#########################################################
#' Computes between consecutive entries of a vector,
#' the differences.
#'
#' @param x {vector} a numeric vector (length n)
#'
#' @return {vector} the differences between consecutive entries (length n-1)
Math.__differencesBetweenNeighbouredValues <- function(x) {
  xPlusi <- x[-1];
  xi <- x[-length(x)];
  
  return(xPlusi - xi);
}


#########################################################
#' Translates signs of values into rewards (G-variable values) 
#' of the Gleichlauf-formula.
#'
#' @param value {numeric} the value which should be converted
#'
#' @return {numeric} the translated value
#' @export
Math.__translateToReward <- function(value) {
  if (value > 0) {
    return(0.5);
  } else if (value < 0) {
    return(-0.5);
  } 
  return(0);
}


#########################################################
#' Computes t-values with a modified t-Test 
#' (T-Test with Hollstein Detrending).
#'
#' @param x {vector} a numeric vector (length(x) > 3)
#' @param y {vector} a numeric vector of same size as x
#'
#' @return {vector} the t-value
#' @source Crossdating_measures.pdf from 
#' Chair of Forest Growth and Dendroecology
#' @export
Math.twoSampleTtestValues <- function(x, y) {
  xDetrended <- Math.__detrendBetweenNeighbouredValues(x);
  yDetrended <- Math.__detrendBetweenNeighbouredValues(y);
  
  correlationCoefficient <- cor(xDetrended, yDetrended);
  tValue <- (correlationCoefficient * sqrt(length(x)-2))/ (1-correlationCoefficient^2);
  
  return(tValue);
}


#########################################################
#' Computes Hollstein Wuchswerte 
#'
#' @param x {vector} a numeric vector (length n)
#'
#' @return {vector} the detrend-values
#' @source Crossdating_measures.pdf from 
#' Chair of Forest Growth and Dendroecology
Math.__detrendBetweenNeighbouredValues <- function(x) {
  xPlusi <- x[-1];
  xi <- x[-length(x)];
  
  return(log(xi / xPlusi));  # natural log
}