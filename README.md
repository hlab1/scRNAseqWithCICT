# CICT for Gene Regulatory Network Inference from Gene Expression Data

**Reference**: Abbas Shojaee and Shao-shan Carol Huang, Robust discovery of gene regulatory networks from single-cell gene expression data using Causal Inference with Composition of Transactions. Manuscript in preparation. 2023.

This repository contains code to run the CICT (Causal Inference Using Composition of Transactions) algorithm on various datasets. The code is organized into two main parts:

* Initial setup and configuration
* Execution of the CICT algorithm

## 1. Initial setup and configuration

Before running the CICT algorithm, you need to set up your working environment, input data files, and CICT configuration file.  We follow the directory structure of the [BEELINE pipeline](https://github.com/Murali-group/Beeline), where `PROJ_ROOT` is the top level directory where this repo is cloned.

### 1a. Install R dependencies

The pre-requisite packages can be installed by running the `requirements.r` script. Or, if you have [renv](https://rstudio.github.io/renv/) set up,  you can call the function `renv::restore()` in R from the `PROJ_ROOT` directory, which will use the `renv.lock` file provided to restore the project library.

### 1b. Make configration file

See example in `config-files/config_SERGIO_DS4_net01_CICT.yaml`.  It defines the paths of the input and output data files and parameters for running CICT.

### 1c. Organize expression data file

The file paths are constructed from values specified in the configuration file in 1b. The gene expression data file is a comma-separated text file in the path  `{input_dir}/{dataset_dir}/{datasets.name}/{exprData}`. The ground-truth network is a comma-seperated text file in `{input_dir}/{dataset_dir}/{datasets.name}/{refNetwork}`.  See examples in `inputs/SERGIO_DS4/net0` and `inputs/SERGIO_DS4/net1`.


## 2. Running CICT

### 2a. Running CICT on the command line

CICT can be called directory on the command line with the configuration file and command line arguments in this format:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R <operation> <config_file_path> <force_output> [<use_preset_learning>]
```
The arguments are:
* \<operation>\: One of the following options: `calcEdges`, `runCICT`, `runSupervised`, `config_par`, `calcEdges_par`, `runCICT_par` 
* <config_file_path>: Path to the configuration file for the CICT algorithm
* <force_output>: Set to 'TRUE' to overwrite existing outputs, 'FALSE' otherwise
* <use_preset_learning> (optional): Set to `TRUE` to use use existing learning set edges in `train.csv`and `test.csv`, default to `FALSE`


Basic workflow of CICT involves two operations: `calcEdges` and `runCICT`.

First use operation `calcEdges` to calculate raw edge weights:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R calcEdges config-files/config_SERGIO_DS4_net01_CICT.yaml TRUE
```
This creates the file `{input_dir}/{dataset_dir}/{datasets.name}/rawEdges.csv`

Then use the operation `runCICT` to conduct CICT training and prediction using learning sets sampled from the ground truth:
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R runCICT config-files/config_SERGIO_DS4_net01_CICT.yaml TRUE
```
This creates the inferred network in output file `{output_dir}/{dataset_dir}/{datasets.name}/CICT/rankedEdges.csv`.

If you want to use existing training and test set files in the `{output_dir}/{dataset_dir}/{datasets.name}/CICT` folder, add the fourth argument (`use_preset_learning`) set it to TRUE. 
```
cd $PROJ_ROOT
Rscript Algorithms/CICT/runCICTEval2.R runCICT config-files/config_SERGIO_DS4_net01_CICT.yaml TRUE TRUE
```
This uses the `train.csv` and `test.csv` learning sets from the previous run.

### 2b. Running CICT in R using configuration file and `runCICTEval2.R` driver script

Start R in the `PROJ_ROOT` directory.

Set the varaible `args.cmnd` for calculating raw edge weights or run CICT (see 2a) and source the driver script `runCICTEval2.R`.

To calculate edge weights:
```
args.cmnd = c('calcEdges','config-files/config_SERGIO_DS4_net01_CICT.yaml', TRUE) 
source('Algorithms/CICT/runCICTEval2.R')
```

To run CICT training and prediction with new sampling of training and test sets:
```
args.cmnd = c('runCICT','config-files/config_SERGIO_DS4_net01_CICT.yaml', TRUE) 
source('Algorithms/CICT/runCICTEval2.R')
```

To run CICT training and prediction with training and test sets from a previous run:
```
args.cmnd = c('runCICT','config-files/config_SERGIO_DS4_net01_CICT.yaml', TRUE, TRUE)
source('Algorithms/CICT/runCICTEval2.R')
```


## 3. Additional Notes

* If you want to run the CICT algorithm in parallel, use the future.batchtools and batchtools libraries and configure the parallel settings as needed.

* You can modify the code to add new datasets, change the edge types, and adjust the algorithm settings as needed.  The code also supports sensitivity analysis and scaling tests for evaluating the performance of the CICT algorithm. See the operations `config_par`, `calcEdges_par`, `runCICT_par`.

* Please refer to the code comments for more information on the CICT algorithm and its configuration options.

## 4. License
scRNAseqWithCICT is licensed under Creative Commons Attribution-NonCommercial-ShareAlike (BY-NC-SA) (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode) provided the following addendum:

1. **Definitions**

   1.1. **Source Code:** The preferred form of the work for making modifications to it.

   1.2. **Object Code:** Any non-source form of a work.

   1.3. **Covered Work:** The work of authorship, whether in Source or Object form, conveyed under this License.

   1.4. **Network Use:** Use of the Covered Work by making it available over a network to the public or specific individuals, whether the network is public or private, regardless of the technology used.

   1.5. **Commercial Use:** Any use of the Covered Work that is primarily intended for or directed towards commercial advantage or monetary compensation. 

2. **Grant of Rights**

   2.1. **Copyright and Related Rights:** As delineated by the original CC BY-NC-SA License.

   2.2. **Patent Rights:** Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims in its contribution, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contribution.

3. **Limitations on Rights**

   3.1. **Noncommercial Use:** You may copy, distribute, display, perform, and modify the Covered Work, but solely for Noncommercial purposes. For any Commercial Use of the Covered Work, you must obtain separate permission from the rights holder.

   3.2. **Network Use:** If you make the Covered Work available over a network, you must also provide the Source Code, so that users interacting with the work over the network may also modify or create derivatives of the Covered Work under the same terms as in this License.

4. **Disclaimer of Warranty**

   4.1. Unless separately stated in writing by the copyright holder, the Covered Work is provided under this License on an "as is" basis, without warranties or conditions of any kind, either express or implied including, without limitation, any warranties or conditions of title, non-infringement, merchantability, or fitness for a particular purpose. 

   4.2. You are solely responsible for determining the appropriateness of using or redistributing the Work and assume any risks associated with Your exercise of permissions under this License.
