#' University of Freiburg SS 2018
#' Chair for Bioinformatics
#' Supervisors: Dr. Martin Raden, PD Dr. Hans-Peter Kahle
#' Authors: Alexander Mattheis

compiler::enableJIT(3);  # compile functions (closures), loops and more

CLUSTERING_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file

#' Allows to optimally cluster data.

#########################################################
#' Creates a small subset of the curves 
#' in which the variance is obtained by measuring
#' the distance of first curve other curves. 
#' Therefore one normalized curve is extracted and from this curve 
#' the distance to all other normalized curves is measured.
#' An amount i.e. 12 of the curves with lowest produced distances 
#' and biggest produced distances are added to the new profilesToAlign 
#' and profilesForAlignment back.
#'
#' @param profilesToAlign {list} the normalized curves which are aligned
#' @param profilesForAlignment {list} the curves on which the alignment is applied to compute a consensus
#' 
#' @return {list(normalizedProfiles, newProfilesForAlignment)} the subsets
Clustering.getSubsetOfProfiles <- function(profilesToAlign, profilesForAlignment) {
  firstCurve <- profilesToAlign[[1]];
  
  # remove first profile
  profilesToAlignMinus1 <- profilesToAlign[-1];
  profilesForAlignmentMinus1 <- profilesForAlignment[-1];
  
  # iterate over all curves and compute mean-distances
  meanDistances <- list();
  
  for (i in 1:length(profilesToAlignMinus1)) {
    nextCurve <- profilesToAlignMinus1[[i]];
    meanDistance <- Math.getMeanDistance(firstCurve, nextCurve);
    meanDistances <- rlist::list.append(meanDistances, meanDistance);
  }
  
  # sort curves by mean distances
  order <- order(unlist(meanDistances));
  
  sortedProfilesToAlign <- profilesToAlignMinus1[order];
  sortedProfilesForAlignment <- profilesForAlignmentMinus1[order];
  
  return(Clustering.getBestAndWorst(Parameters.SELECTION_NUMBER, sortedProfilesToAlign, sortedProfilesForAlignment));
}


#########################################################
#' Returns a number of non-similar and similar curves.
#'
#' @param number {number} the number of non-similar and similar profiles to select
#' @param sortedProfilesToAlign {list} the sorted normalized curves which are aligned
#' @param sortedProfilesForAlignment {list} the sorted curves on which the alignment is applied to compute a consensus
#'
#' @return {list(normalizedProfiles, newProfilesForAlignment)} the subsets
Clustering.getBestAndWorst <- function(number, sortedProfilesToAlign, sortedProfilesForAlignment) {
  # return best and worst
  addedGoodInd <- 1:Parameters.SELECTION_NUMBER;
  addedBadInd <- (length(sortedProfilesToAlign) - number + 1):length(sortedProfilesToAlign);
  
  newProfilesToAlign <- append(sortedProfilesToAlign[addedGoodInd], sortedProfilesToAlign[addedBadInd]);
  newProfilesForAlignment <- append(sortedProfilesForAlignment[addedGoodInd], sortedProfilesForAlignment[addedBadInd]);
  
  return(list(normalizedProfiles = newProfilesToAlign, profiles = newProfilesForAlignment));
}


#########################################################
#' Creates a smaller subset of the curves in which
#' the variance is obtained by measuring
#' the all-vs-all average distances between curves.
#' Therefore a normalized curve is extracted and from this curve 
#' the distance to all other normalized curves is measured.
#' An amount i.e. 12 of the curves with lowest produced distances 
#' and biggest produced distances are added to the new profilesToAlign 
#' and profilesForAlignment back. 
#' So most similar and most different are added together.
#' Hint: The computation was very fast such that the loading 
#' from the harddrive was omitted.
#' 
#' @param profilesToAlign {list} the normalized curves which are aligned
#' @param profilesForAlignment {list} the curves on which the alignment is applied to compute a consensus
#' 
#' @return {list(normalizedProfiles, newProfilesForAlignment)} the subsets
Clustering.getSubsetOfProfiles2 <- function(profilesToAlign, profilesForAlignment) {
  order <- Clustering.getOrder(profilesToAlign);
  
  sortedProfilesToAlign <- profilesToAlign[order];
  sortedProfilesForAlignment <- profilesForAlignment[order];
  
  return(Clustering.getBestAndWorst(Parameters.SELECTION_NUMBER, 
                                       sortedProfilesToAlign, 
                                       sortedProfilesForAlignment));
}


########################################################
#' Returns the order for profiles depending on their mean-distance.
#'
#' @param profiles {list} the profiles order dependent on mean distance
#'
#' @return {vector} the order depending on the mean-distances
#' @export
Clustering.getOrder <- function(profiles) {
  # create empty distance matrix initialized with diag = 0
  numProfiles <- length(profiles);
  mat <- matrix(nrow = numProfiles, ncol = numProfiles);
  diag(mat) <- 0;
  
  # iterate over all profiles in that year and compute pairwise distances
  for (i in 1:length(profiles)) {
    profile <- profiles[[i]];
    
    distances <- Analyzer.getPairwiseDistances(profiles[-(1:i)], profile, TRUE);
    
    if (length(distances) > 0) {
      # iterate over all inner list
      for (j in 1:length(distances)) {
        value <- distances[[j]];
        
        mat[[i, j+i]] <- value;  # +1 since first position is zero
        mat[[j+i, i]] <- value;
      }
    }
  }
  
  # compute meanDistances for each column of the matrix
  numMeanDistances <- numProfiles - 1;  # n profiles -> n-1 distances
  meanDistances <- colSums(mat) / numMeanDistances;
  
  # sort profiles with that distances
  order <- order(unlist(meanDistances));
  
  return(order);
}


#########################################################
#' Clusters the profiles in each year 
#' with the gap statistic approach and k-Means. 
#' Then it selects for each year an amount of profiles 
#' by uniform selection of a single profile from each cluster.
#' Hint: It was too unstable!
#'
#' @param perYearProfiles {enviroment} the per year profiles which should be clustered
#' @param startYear {numerical} the year at which the chronology starts
#' @param endYear {numerical} the year at which the chronology ends
#' @param limit {numerical} the minimum number of profiles necessary to start a clustering
#' @param selectionNumber {numerical} the number of profiles to select for final bucket
#' (for value 1, it is selected the 1 profile from each cluster)
#' 
#' @return {environment} the environment which contains maximum limit many profiles per year
Clustering.getPerYearProfilesPerCluster <- function(perYearProfiles, startYear, endYear, limit, selectionNumber) {
  profilesPerYear <- new.env();
  
  print(Titles.CURRENTLY_PROCESSED_YEAR);
  
  # iterate over each year
  for (i in startYear:endYear) {
    print(i);
    profiles <- perYearProfiles[[toString(i)]];
    
    clustersProfiles <- list();
    
    # limit the number of profiles
    if (length(profiles) >= limit) {
      clustersProfiles <- Clustering.__getProfilesFromGapStatisticClustering(profiles, numBootstraps = 30, 
                                                                             selectionNumber = selectionNumber);
    } else {
      clustersProfiles <- profiles;
    }
    
    profilesPerYear[[toString(i)]] <- clustersProfiles;
  }
  
  return(profilesPerYear);
} 


#########################################################
#' After choosing the optimal number of clusters 
#' for k-means, from each cluster the profile 
#' which is the most similar to the others returned.
#'
#' @param profiles {list} the profiles which have to be clustered
#' @param numBootstraps {numerical} the number of bootstraped distributions
#' @param selectionNumber {numerical} the number of profiles to select for final bucket
#' (for value 1, it is selected the 1 profile from each cluster)
#' 
#' @return {list} the one profile per cluster
#' @export
Clustering.__getProfilesFromGapStatisticClustering <- function(profiles, numBootstraps, selectionNumber) {
  X <- matrix(unlist(profiles), nrow = length(profiles), byrow = TRUE);
  clustering <- Clustering.__getOptimalKvalueClustering(X, numBootstraps);

  # create clusters containing profiles
  clusters <- list();
  uniqueClusters <- unique(clustering);
  
  # print(length(uniqueClusters));
  
  for (i in 1:length(uniqueClusters)) {
    clusterName <- uniqueClusters[[i]];
    indicesInX <- which(clustering == clusterName);
    
    cluster <- list();
    
    for (j in indicesInX) {
      profile <- X[j, ];
      cluster <- rlist::list.append(cluster, profile);
    }
    
    clusters <- rlist::list.append(clusters, cluster);
  }
  
  profiles <- list();
  
  stop <- FALSE;
  
  if (selectionNumber == 1) {
    # selects the profile with the lowest distance to the others in the cluster
    for (i in 1:length(clusters)) {
      cluster <- clusters[[i]];
      bestIndex <- Analyzer.__getIndexOfMostSimilarProfile(cluster);
      profile <- as.list(cluster[[bestIndex]]);
      profiles <- rlist::list.append(profiles, profile);
    }
  } else {  # get selection number many profiles
    while (length(profiles) < selectionNumber) {
      if (stop) { break; }
      
      clusterIndices <- sample(1:length(clusters), length(clusters));  # to get indices in a random order
      
      for (i in clusterIndices) {
        cluster <- clusters[[i]];
        
        if (length(cluster) > 0) {
          index <- sample(1:length(cluster), 1);
          profile <- as.list(cluster[[index]]);
          profiles <- rlist::list.append(profiles, profile);
          cluster <- cluster[-index];
          clusters[[i]] <- cluster;
        }
        
        if (length(profiles) >= selectionNumber) {
          stop <- TRUE;
          break;
        }
      }
    }
  }

  return(profiles);
}


#########################################################
#' Chooses best value k with gap statistic test for k-Means and returns the corresponding clustering.
#'
#' @param X {matrix} the matrix with profile points as features
#' @param numBootstraps {numerical} the number of bootstraped distributions
#'
#' @return {vector} a vector telling to which cluster the different profiles were added to
Clustering.__getOptimalKvalueClustering <- function(X, numBootstraps) {
  clusterings <- list();
  dispersions <- list();
  numClustersMax <- nrow(X)-1;
  
  # [1] execute k-means with different number of clusters
  for (i in 1:numClustersMax) {
    clustering <- kmeans(X, i, nstart = 10, iter.max = 250);
    dispersion <- clustering$tot.withinss;
    
    clusterings <- rlist::list.append(clusterings, clustering);
    dispersions <- rlist::list.append(dispersions, dispersion);
  }
  
  # [2] compute B reference sets (with same amount of points) and gap statistic
  bootstrapDispersions <- list();
  
  for (b in 1:numBootstraps) {
    features <- list();
    
    # iterate over each column and create random feature-points of same length
    for (j in 1:ncol(X)) {
      # get range for column
      column <- X[, j];
      
      minimum <- min(column);
      maximum <- max(column);
      
      newColumn <- runif(nrow(X), minimum, maximum);
      features <- rlist::list.append(features, newColumn);
    }
    
    newX <- do.call(cbind, features);
    
    dispersionsPerClusterSize <- list();
    
    for (i in 1:numClustersMax) {
      clustering <- kmeans(newX, i, nstart = 10, iter.max = 250);
      dispersion <- clustering$tot.withinss;
      dispersionsPerClusterSize <- rlist::list.append(dispersionsPerClusterSize, dispersion);
    }
    
    bootstrapDispersions <- rlist::list.append(bootstrapDispersions, dispersionsPerClusterSize);
  }
  
  # compute gap statistic values
  gapValues <- list();
  averageValues <- list();
  
  for (k in 1:numClustersMax) {
    # get bootstraped average value
    sum <- 0;
    
    for (b in 1:numBootstraps) {
      bootstrapDispersion <- bootstrapDispersions[[b]][[k]];
      sum <- sum + log(bootstrapDispersion);
    }
    
    averageValue <- sum/numBootstraps;
    
    # get final value
    gapValue <- averageValue - log(dispersions[[k]]);
    gapValues <- rlist::list.append(gapValues, gapValue);
    averageValues <- rlist::list.append(averageValues, averageValue);
  }
  
  # [3] compute standard deviations and simulation errors
  simulationErrors <- list();
  
  for (k in 1:numClustersMax) {
    sum <- 0;
    
    for (b in 1:numBootstraps) {
      bootstrapDispersion <- bootstrapDispersions[[b]][[k]];
      sum <- sum + (log(bootstrapDispersion) - averageValues[[k]])^2;
    }
    
    sd <- sqrt(sum/numBootstraps);
    simulationError <- sd * sqrt(1 + 1 / numBootstraps);
    simulationErrors <- rlist::list.append(simulationErrors, simulationError);
  }
  
  # [4] choose smallest k, which fullfills gap(k) >= gap(k+1) - s_{k+1}, else set k = 1
  smallestK <- 1;
  
  for (k in 2:(numClustersMax-1)) {
    if (gapValues[[k]] >= gapValues[[k+1]] - simulationErrors[[k+1]]) {
      smallestK <- k;
      break;
    }
  }
  
  # View(gapValues);
  # print(smallestK);
  
  return(clusterings[[smallestK]]$cluster);
} 