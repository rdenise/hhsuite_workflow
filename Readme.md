# Snakemake workflow: HHsuite workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.14.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/rdenise/hhsuite_workflow/workflows/Tests/badge.svg?branch=main)](https://github.com/rdenise/hhsuite_workflow/actions?query=branch%3Amain+workflow%3ATests)
[![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](LICENSE)

## Aim

This software will produce a plot of the similarities between the hmm profiles created from the alignment files that you give as input. The workflow will produce multiples files using hhsuite package. The first one is the hmm profiles in HHM format using hhmake, the second is the results of hhblits that will compare the hmm profiles between each other. (results will be in hhr and tsv format). It will also produce the hmm profiles in HMM format using hhbuild from hmmer package.

## Installation

### Step 1: install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 

#### Notes 

For all following commands (step 2 to 5) ensure that this environment is activated via 

```shell
conda activate snakemake
```

### Step 2: deploy workflow

 Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system in the place of your choice as follow (note to change the path and file name to the one you want to create): : 

```shell
mkdir path/to/project-workdir
```

Then go to your project working directory as follow:

```shell
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.

Second, run 

```shell
snakedeploy deploy-workflow https://github.com/rdenise/hhsuite_workflow . --tag 0.1.3
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

### Step 3: configure workflow

#### General settings

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.  

annotations sheet
- (Optional) if you want to define your own color, add alignment name to `config/annotations.tsv`. For each protein, the columns `HMM`, and `color` have to be defined. This file is optional, you can put the value to '' if you don't have one and want the workflow to create the color for you. The `HMM` column correspond to the name of the alignment file. The `color` columns correspond to the color you want in HEX format. Example of the `annotations file` is present in the `config` folder if needed or in the [doc](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_annotations.tsv) folder in the GitHub page

### Step 4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda, run Snakemake with 

```shell
snakemake --cores 1 --use-conda 
```

### Step 5: generate report

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via 

```shell
snakemake --report report.zip --report-stylesheet config/report.css
```
The resulting report.zip file can be passed on to collaborators, provided as a supplementary file in publications.


## Walk-Through and File Production

This pipeline consists of 6 steps called rules that take input files and create output files. Here is a description of the pipeline.

1. As snakemake is set up, there is a last rule, called `all`, that serves to call the last output files and make sure they were created.

2. A folder containing your work will be created:

```
   [project_name]/                           <- top-level project folder (your project_name)
   │
   │
   ├── report                                <- Folder with the report.html file inside    
   │
   ├── logs                                  <- Collection of log outputs
   │
   ├── databases                             <- Generated analysis database related files
   │    
   ├── alignments                            <- If conversion needed the a3m format of the alignment will be here
   │
   ├── hhblits                               <- Folder with the results of hhblits
   │   ├── hhr                               <- Results of hhblits in hhr format
   │   └── tsv                               <- Results of hhblits in tsv format
   │
   ├── hhm                                   <- HMM profiles in HHM format made with hhmake
   ├── hmm                                   <- HMM profiles in HMM format made with hmmbuild
   │
   └── plots                                 <- Folder of the final results with the figure in pdf and png 


```

### Pipeline in image 

#### Normal behavior

<p align="center">
  <img src="doc/dummy_dag.png?raw=true" height="400">
</p>

