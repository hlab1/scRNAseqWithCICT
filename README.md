
# CICT for Gene Regulator Network Infernce from gene expression data

This repository contains code to run the CICT (Causal Inference Causality Test) algorithm on various datasets. The code is organized into two main parts:

    Initial setup and configuration
    Execution of the CICT algorithm

Getting Started

Before running the CICT algorithm, you need to set up your working environment and specify the paths for input and output files.
Step 1: Set working directory and paths

    Set the working directory by modifying the url.base variable to point to the base directory where the input files and folders are located.
    Set the paths for input files, ground truth files, and output folders by modifying the corresponding variables in the code.

Step 2: Load required libraries and functions

    Source the required libraries and functions by running the following line:

R

source('/scratch/as15096/eric/Algorithms/CICT/requirements/CICT_LibsFunctions.R')

Step 3: Configure the CICT algorithm

    The args.cmnd variable contains the command line arguments for the CICT algorithm. It should be in the following format:

R

args.cmnd = c('operation', 'config_file_path', 'force_output')

where:

    operation: One of the following options: 'calcEdges', 'runCICT', 'runSupervised', 'runCICT_par', 'install', 'calcEdges_par'
    config_file_path: Path to the configuration file for the CICT algorithm
    force_output: Set to 'TRUE' to overwrite existing outputs, 'FALSE' otherwise

Running the CICT Algorithm

After setting up your environment and specifying the necessary paths and configurations, you can run the CICT algorithm using the following steps:

    Run the CICT function by providing the required input parameters:

R

CICT(theJobID, url.input, url.rawedgefile, url.name.map, url.gt, url.output, url.logfile, cictRawEdgeCol, earlyThresholdForGraphAnalysis, minGroundTruth.ratio.learning, maxunseenTest.ratio, maxGroundTruth, randomEdgesFoldCausal, sampling.c_rc.ratio, trainingTarget, tstPrecent, forceOutput, arg.experiment, FLAG_runOnAllEdges, FLAG_exportRankedEdges, FLAG_exportTrainAndTest, Debug, RF_max_depth, RF_ntrees, preset.train, preset.test, maxNetSize, ...)

    If you want to run the CICT algorithm in parallel, use the future.batchtools and batchtools libraries and configure the parallel settings as needed.

Example Configuration

An example configuration for running the CICT algorithm is provided below:

R

args.cmnd = c('runCICT_par', '/scratch/as15096/eric/outputs/cict_par/sens_sparsity/parConf_7.yaml', 'FALSE')

This configuration will run the CICT algorithm in parallel using the specified configuration file and will not overwrite existing outputs.
Additional Notes

    You can modify the code to add new datasets, change the edge types, and adjust the algorithm settings as needed.
    The code also supports sensitivity analysis and scaling tests for evaluating the performance of the CICT algorithm.

Please refer to the code comments for more information on the CICT algorithm and its configuration options.
