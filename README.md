
# CICT for Gene Regulator Network Infernce from gene expression data

This repository contains code to run the CICT (Causal Inference Causality Test) algorithm on various datasets. The code is organized into two main parts:

* Initial setup and configuration
* Execution of the CICT algorithm

## 1. Initial setup and configuration

Before running the CICT algorithm, you need to set up your working environment, input data files, and CICT configuration file.  We follow the directory structure of the [BEELINE pipeline](https://github.com/Murali-group/Beeline), where `PROJ_ROOT` is the top level directory where this repo is cloned.

### 1a. Install R dependencies

The pre-requisite packages are installed by the `requirement.r` script. 

### 1b. Make configration file

See example in `config-files-split/config_L0_split/hHep/CICT/config.yaml`.  It defines the paths of the input and output data files and parameters for running CICT.

### 1c. Organize expression data file

The file paths are constructed from values specified in the configuration file in 1b. The gene expression data file is a comma-separated text file with header in the path  `{input_dir}/{dataset_dir}/{datasets.name}/{exprData}`. The ground-truth network is a comma-seperated text file with header in `{input_dir}/{dataset_dir}/{datasets.name}/{refNetwork}`.  See examples in `inputs/L0/hHep/ExpressionData.csv` and `inputs/L0/hHep/refNetwork.csv`.


## 2. Running CICT

### 2a. Running CICT on the command line

CICT can be called directory on the command line with the configuration file with command line arguments in this format:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R <operation> <config_file_path> <force_output> [<use_preset_learning>]
```

\<operation\>: One of the following options: `calcEdges`, `runCICT`, `runSupervised`, `config_par`, `calcEdges_par`, `runCICT_par` 

<config_file_path>: Path to the configuration file for the CICT algorithm
  
<force_output>: Set to 'TRUE' to overwrite existing outputs, 'FALSE' otherwise
  
<use_preset_learning> (optional): Set to `TRUE` to use use existing learning set edges in `train.csv`and `test.csv`, default to `FALSE`


* To run CICT, first use operation `calcEdges` to calculate raw edge weights:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R calcEdges config-files-split/config_L0_split/hHep/CICT/config.yaml TRUE
```
This creates the file `{input_dir}/{dataset_dir}/{datasets.name}/rawEdges.csv`

* use the operation `runCICT` to conduct CICT training and prediction using learning sets sampled from the ground truth:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R runCICT config-files-split/config_L0_split/hHep/CICT/config.yaml TRUE
```
* Alternatively, to use existing training and test set files in the `{output_dir}/{dataset_dir}/{datasets.name}/CICT` folder, add the `use_preset_learning` argument and set it to TRUE.
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R runCICT config-files-split/config_L0_split/hHep/CICT/config.yaml TRUE TRUE
```
This creates the inferred network in output file `{output_dir}/{dataset_dir}/{datasets.name}/CICT/rankedEdges.csv`.

## 2b. Running CICT in R using configuration file and `runCICTEval2.R` driver script

* Start R in the `PROJ_ROOT` directory.

* Source the `runCICTEval2.R` script to load the required libraries and functions
```
source('Algorithms/CICT/runCICTEval2.R')
```

* Set the varaible `args.cmnd` for calculating raw edge weights or run CICT (see 2a) and call the driver function `runCICTEval2`
```
args.cmnd = c('calcEdges','config-files-split/config_L0_split/hHep/CICT/config.yaml', TRUE) 
runCICTEval2(args.cmnd)
```


### 2c. Running CICT in R using `CICT` fucntion
The `CICT` function allows the CICT to be run without a configuration file.
```
CICT(theJobID, url.input, url.rawedgefile, url.name.map, url.gt, url.output, url.logfile, cictRawEdgeCol, earlyThresholdForGraphAnalysis, minGroundTruth.ratio.learning, maxunseenTest.ratio, maxGroundTruth, randomEdgesFoldCausal, sampling.c_rc.ratio, trainingTarget, tstPrecent, forceOutput, arg.experiment, FLAG_runOnAllEdges, FLAG_exportRankedEdges, FLAG_exportTrainAndTest, Debug, RF_max_depth, RF_ntrees, preset.train, preset.test, maxNetSize, ...)
```

## 3. Additional Notes

* If you want to run the CICT algorithm in parallel, use the future.batchtools and batchtools libraries and configure the parallel settings as needed.

* You can modify the code to add new datasets, change the edge types, and adjust the algorithm settings as needed.  The code also supports sensitivity analysis and scaling tests for evaluating the performance of the CICT algorithm. See the operations `config_par`, `calcEdges_par`, `runCICT_par`.

* Please refer to the code comments for more information on the CICT algorithm and its configuration options.
