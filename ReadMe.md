<img src="https://github.com/AlexanderMattheis/Cross-Dating/blob/master/crossdating_icon_large.png" width="128" height="128">

# Cross-Dating of Intra-Annual<br/>Wood Density Series
<a name="introduction"></a>
The age of wood pieces can be determined by annual rings. 
You can just measure the widths or maximum densities of individual rings 
and then shift them in a chronology with known dates to determine the exact age. 
But such established approaches do not work well with short pieces, 
i.e. with wood pieces containing only a few annual rings. Our new approach, however, is different. 
It has been shown that it can accurately determine the date for even shorter wood pieces correctly. 
Therefore, it uses the series of densities within a ring.
By this it becomes one of the most accurate, existing approaches in dendrochronology, today!

## Overview
- [Introduction](#introduction)
- [Installation](#installation)
	- [Quick-Test](#quick-test)
- [Interface](#interface)
	- [Input](#input)
	- [Output](#output)
	- [Provided Functions](#provided-functions)
- [Conventions](#conventions)
	
## Installation
It has been used different libraries for visualization, simplification, testing and fast code execution.

For visualizations and unit-testing following libraries are necessary:

| Library     					                               | Version   | Description
|:-------------------------------------------------------------|:----------|:-----------------------------------------------|
| [ggplot2](http://ggplot2.org/)                               | 2.2.1     | grammar of graphics based plotting system      |
| [gridExtra](https://github.com/cran/gridExtra)               | 2.3 	   | arranging ggplots in grids                     |
| [reshape2](https://github.com/hadley/reshape)                | 1.4.3     | fast reshaping with C++ calls                  |
| [scales ](https://github.com/hadley/scales)                  | 0.4.1     | custom axes for ggplot                         |
| [testthat](http://testthat.r-lib.org/)                       | 2.0.0     | unit-testing maths, loading and visualization  |
| [VennDiagram](https://doi.org/10.1186/1471-2105-12-35)       | 1.6.2     | to create Venn- and Euler-diagrams             |

These libraries can be installed with the command 
`install.packages(c("ggplot2", "gridExtra", "reshape2", "scales", "testthat", "VennDiagram"))`,
but also each package in the table is linked to its repository, paper or homepage.
Perhaps you have to install the right version for correct execution of the code!

Essential libraries for fast and nice code:

| Library                                                   |Version    | Description                                                       |
|:----------------------------------------------------------|:----------|:------------------------------------------------------------------|
| [fitdistrplus](https://github.com/cran/fitdistrplus)      | 1.09      | fitting distributions to data	            						|
| [plyr](https://github.com/hadley/plyr)                    | 1.8.4     | data reconversion with C++ calls    								|
| [rlist](https://github.com/renkun-ken/rlist)              | 0.4.6.1   | to build complex data-structures and avoid unnecessary code       |
| [stringr](https://github.com/tidyverse/stringr)           | 1.2.0     | for consistent usage of string-operations  						|


These libraries can be installed with the command 
`install.packages(c("fitdistrplus", "plyr", "rlist", "stringr"))`.
Perhaps you have to install the right version for correct execution of the code!

In **addition**, MICA has to be installed.
A manual, therefore, can be found [here](https://github.com/BackofenLab/MICA#installation).
It should be enough to execute the command `install.packages("rJava")`.
The current Version of MICA in the `libraries`-folder is 2.02.
Therefore, you have install the Java runtime.
[Java 8 Update 152](https://www.java.com/) together with [R 3.4.3](https://www.r-project.org/)
were successfully tested within this project.

### Quick-Test
The different approaches can be easily tested
after you have installed the necessary libraries.
Set up the working directory to the folder containing
the file called `Main.R`.
The function `Main.interface()` which is executed
by the `Main.main`-function contains preset examples to test the approaches.
Uncomment with `Ctrl+Shift+C` there the selected lines of an approach you want execute. 
Then select all code with `Ctrl+A` in the `Main.R` and press `Ctrl+Enter`.
The execution takes a few minutes. 
Meanwhile, the finished number of samples and 
different information is written into the console to give the user a visual feedback.

To execute the unit-tests, just go into the folder `tests`
and execute the file `run_tests.R`. Do **not** set the folder as the working directory!
	
## Interface
The interface allows you executing each 
base approach presented in the theoretical part.
Hereby you can date samples of inter-annual wood density-series.
By this you get as an output a matrix with most probable dates 
together with corresponding reached scores and p-values.

### Input
The inputs are `*.csv`-files with following structure:

| year   | density | characteristic  |
|:------:|:-------:|:---------------:|
| 1992   | 2.016   | 166             |
| 1992   | 2.433   | 166             |
| 1992   | 2.881   | 166             |
| 1993   | 2.043   | 128             |
| 1993   | 2.383   | 128             |
| &#8942;| &#8942; | &#8942;         |

The rows with the same year correspond to a density profile.
Optionally available, depending on the approach is the `characteristic`-column.
It is implemented an approach with two-steps which can use ring-widths or maximum densities
to speed up the dating procedure. Therefore, it uses only the
most appropriate years from the ring-width or maximum densities based approach.
And there this `characteristic`-column is used to store ring-widths or maximum densities.  
Samples you want to date, must have a similar format.
Concretely the column `year` is replaced by a column called `part`
with numbers 1,2,3,... identifying the different profiles.
For testing purposes, a chronology and samples in the presented formats 
were prepared under `\input\interface\pass_1`.

### Output
The output you get by interface functions is a `matrix`
with the following structure:

| sample 			| pValue      | rank1   | score1  | rank2   | &#8230; |
|:-----------------:|:-----------:|:-------:|:-------:|:-------:|:-------:|
| 1041_MICA-cons    | 0.06209708  | 1960    | 41.6288 | 1947    | &#8230; |
| 1051_MICA-cons    | 0.01127934  | 1987    | 14.9864 | 1940    | &#8230; |
| 1201_MICA-cons_1  | 0.06737666  | 1957    | 24.6861 | 1955    | &#8230; |
| 1201_MICA-cons_2  | 0.066093    | 1941    | 24.6395 | 1962    | &#8230; |
| &#8942;           | &#8942;     | &#8942; | &#8942; | &#8942; | &#8230; |

whereas the `pValue`-column as well as the `score`-columns are 
only optionally available in the approaches Two-Step and Bucket.
In the first column you can see the name of the dated sample, beneath
you see the p-value for the rank 1 rated samples.
Then a column with the rank 1 prediction and a column with the corresponding
score. The user can specify the number of predictions, so there can also be
rank 3 or rank 4 predictions.
This matrix can automatically be stored as `*.csv`-file,
so there is an option in the `Interface`-class to store computed data in the `output`-folder.

### Provided Functions
The project provides a file called `Main.R`. It contains a function
`Main.interface` and this function contains examples for each of
the four approaches.

- [`Interface.computeDatesConsensusApproach(..)`](#Interface.computeDatesConsensusApproach): executes the Consensus Approach per sample
- [`Interface.computeDatesBucketApproach(..)`](#Interface.computeDatesBucketApproach): executes the Bucket Approach per sample
- [`Interface.computeDatesPerTreeApproach(..)`](#Interface.computeDatesPerTreeApproach): executes the Per-Tree Approach per sample
- [`Interface.computeDatesVotingApproach(..)`](#Interface.computeDatesVotingApproach): executes the Voting Approach per sample
- [`Interface.computeDatesTwoStepApproach(..)`](#Interface.computeDatesTwoStepApproach): executes first a fast Points-Based Approach (correlation coefficient / t-value based approach) and afterwards the Bucket Approach on the number of potentially correct years found by the Points-Based Approach<br/>

**Hint:**	Computed chronologies and samples must not contain gaps. Check that for each year at least one consensus-profile or bucket is available!

-------------------
<a name="Interface.computeDatesConsensusApproach"></a>
#### `Interface.computeDatesConsensusApproach(consensusPath, consensusName, samplesPath, scoreType, bestYearsMax, save, fileName)`

`Interface.computeDatesConsensusApproach(..)` executes the Consensus Approach per sample.

*Input parameters:* 
- `consensusPath {string}` : the path to the consensus
- `consensusName {string}` : the name of the consensus file
- `samplesPath {string}` : the path to the samples which should be dated
- `scoreType {string}` : the score-type which should be used for computation (`"a"` = y-based, `"b"` = slope-based, `"c"` = z-scores y-based, `"d"` = z-scores slope-based)
- `bestYearsMax {numeric}` : tells how many best ranked years should be stored
- `save {logical}` : tells if the list of possible dates should be stored in the project `output`-folder
- `fileName {string}` : the filename without extension for the stored per sample dates-file

*Output:* 
- `{matrix}` the matrix of possible dates (see [Output](#output))

-------------------
<a name="Interface.computeDatesPerTreeApproach"></a>
#### `Interface.computeDatesPerTreeApproach(consensusPath, consensusName, samplesPath, scoreType, bestYearsMax, save, fileName)`

`Interface.computeDatesPerTreeApproach(..)` executes the Per-Tree Approach per sample.

*Input parameters:* 
- `consensusPath {string}` : the path to the consensus
- `consensusName {string}` : the name of the consensus file
- `samplesPath {string}` : the path to the samples which should be dated
- `scoreType {string}` : the score-type which should be used for computation (`"a"` = y-based, `"b"` = slope-based, `"c"` = z-scores y-based, `"d"` = z-scores slope-based)
- `bestYearsMax {numeric}` : tells how many best ranked years should be stored
- `save {logical}` : tells if the list of possible dates should be stored in the project `output`-folder
- `fileName {string}` : the filename without extension for the stored per sample dates-file

*Output:* 
- `{matrix}` the matrix of possible dates (see [Output](#output))

-------------------
<a name="Interface.computeDatesBucketApproach"></a>
#### `Interface.computeDatesBucketApproach(bucketsPath, samplesPath, scoreType, innerFunc, outerFunc, bestYearsMax, qualityMeasures = "", save, fileName))`
`Interface.computeDatesBucketApproach(..)` executes the Bucket Approach per sample.<br/> 

*Hint:* The chronology must have at least length 30 to get proper p-values.

*Input parameters:* 
- `bucketsPath {string}` : the path to the per tree consensi which should be used
- `samplesPath {string}` : the path to the samples which should be dated
- `scoreType {string}` : the score-type which should be used for computation (`"a"` = y-based, `"b"` = slope-based, `"c"` = z-scores y-based, `"d"` = z-scores slope-based)
- `innerFunc {function}` : the function which should be applied on the per bucket computed scores to get a score for the given position (e.g. the minimum function `min`)
- `outerFunc {function}` : the function which should be applied on the per sample computed scores to get a final score for the position (e.g. the summation function `sum`)
- `bestYearsMax {numeric}` : tells how many best ranked years should be stored
- `qualityMeasures {string}` : tells which quality measures should be active (combine multiple options: `""` = none, `"p"` = p-values, `"s"` = scores, `"ps"` or `"sp"` = scores and p-values)
- `save {logical}` : tells if the list of possible dates should be stored in the project `output`-folder
- `fileName {string}` : the filename without extension for the stored per sample dates-file

*Output:* 
- `{matrix}` the matrix of possible dates

-------------------
<a name="Interface.computeDatesVotingApproach"></a>
#### `Interface.computeDatesVotingApproach(bucketsPath, samplesPath, scoreType, topYearsCount, approach = "", minimumLength, bestYearsMax, save, fileName)`
`Interface.computeDatesVotingApproach(..)` executes the Voting Approach per sample.<br/> 

*Hint 1:* It is possible that `NA`'s are returned for rank predictions, 
since for example all votes could have gone to a single year.<br/>
*Hint 2:* If `minimumLength` is not set correctly the approach will need exponential time.

*Input parameters:* 
- `bucketsPath {string}` : the path to the per tree consensi which should be used
- `samplesPath {string}` : the path to the samples which should be dated
- `scoreType {string}` : the score-type which should be used for computation (`"a"` = y-based, `"b"` = slope-based, `"c"` = z-scores y-based, `"d"` = z-scores slope-based)
- `topYearsCount {numeric}` : the number of years selected per column
- `approach {string}` : the string telling you which approaches should be active (combine multiple options: `""` = none, `"p"` = powerset approach)
- `minimumLength {numeric}` : tells which minimum sample lengths should be considered in the powerset table (`-1` = no limit)
- `bestYearsMax {numeric}` : tells how many best ranked years should be stored
- `save {logical}` : tells if the list of possible dates should be stored in the project `output`-folder
- `fileName {string}` : the filename without extension for the stored per sample dates-file

*Output:* 
- `{matrix}` the matrix of possible dates

-------------------
<a name="Interface.computeDatesTwoStepApproach"></a>
#### `Interface.computeDatesTwoStepApproach(bucketsPath, samplesPath, scoreTypeRingWidths, scoreTypeBuckets, topYearsCount, bestYearsMax, qualityMeasures, save, fileName)`
`Interface.computeDatesTwoStepApproach(..)` executes first a fast Points-Based Approach (correlation coefficient / t-value based approach)
and afterwards the Bucket Approach on the amount of potentially correct years found by the Points-Based Approach.<br/>

*Hint:* `topYearsCount` has to be set at least to 30 to get a proper distribution for p-values!

*Input parameters:* 
- `bucketsPath {string}` : the path to the per tree consensi which should be used
- `samplesPath {string}` : the path to the samples which should be dated
- `scoreTypeCharacteristic {string}` : the score-type which should be used for computation during the Points-Based Approach (`"p"` = Pearson's Rho, `"t"` = Kendall's Tau, `"r"` =  Spearman's Rho, `"v"` = t-value) 
- `scoreTypeBucket {string}` : the score-type which should be used for computation (`"a"` = y-based, `"b"` = slope-based, `"c"` = z-scores y-based, `"d"` = z-scores slope-based)
- `topYearsCount {numeric}` : the number of top years which are stored by the characteristic approach (at least 2)
- `bestYearsMax {numeric}` : tells how many best ranked years should be stored
- `qualityMeasures {string}` : tells which quality measures should be active (combine multiple options: `""` = none, `"p"` = p-values, `"s"` = scores, `"ps"` or `"sp"` = scores and p-values)
- `save {logical}` : tells if the list of possible dates should be stored in the project `output`-folder
- `fileName {string}` : the filename without extension for the stored per sample dates-file

*Output:* 
- `{matrix}` the matrix of possible dates

## Conventions
#### Assignments
R allows the `=`-operator as the assignment-operator, but we use `<-` for assignments instead.


#### Classes
The class name is also the name of the file, and it is written behind the property or function e.g.

```r
Plotter.__extendedMode <- TRUE;

Plotter.getLogNormalDistributionPlot <- function(scores, fit, color, print) {
	...
}
```

like in C++.

#### Constants
Almost all constant values
like paths, strings, symbols, filenames are stored in the `Defaults.R`.
So everything except functional strings 
like `dotted` or `histogram` is stored there.

#### Debugging
To avoid sourcing the same class two times
and by this deactivating breakpoints, there is 
in every class an import-boolean written 
in big letters and named like the class:

```r
EXERCISE_1_IMPORTED <- TRUE;  # to avoid a reimport by the "Main.R"-class after sourcing this file
```
This import-boolean is set, 
when the class is sourced for example
after setting a breakpoint.
If you now execute the `Main`-*class* source code, 
it is automatically checked if the boolean already exists
and a set breakpoint in *class* `Exercise1` won't be deactivated. 
```Exercise1``` is not resourced 
because of an existence check-up in `Main`:

```r
if(!exists("EXERCISE_1_IMPORTED")) source("Exercise1.R");
```

By not using this technique you have otherwise to uncomment 
the sourcing `source("Exercise1.R")` of the *class* `Exercise1` 
in *class* `Main`. And that every time you 
want to debug the *class* `Exercise1`. <br/>

**Hint:**	It is not allowed to have two classes with the same name!


#### Visibility
Visibility follows the Python programming style 
by just marking functions as protected or private. <br/>
Functions with two underscores e.g.

```r
Exercise1.__getSubpatterns <- function(patternY, subpatternsIntervals) {
	...
}
```
are private functions. <br/>

And functions with one underscore e.g.

```r
Alignment._createAlignment <- function(path, sequenceA, sequenceB) {
	...
}
```

are protected functions.