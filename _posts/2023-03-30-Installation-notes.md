---
layout: post
title: "Installation notes"
categories: misc
---

### 1. Install Linux (Ubuntu) using WSL on windows

### 2. Install Miniconda

* First, download the installation [script](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

* Run the installation script:

  ````sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ````
  
* Type ‘yes’ to all questions.

* Close and reopen the shell.

### 3. Install seq-seq-pan

* First, add the correct conda channel to find the software:

  ````sh
  conda config --add channels defaults
  conda config --add channels bioconda
  ````

* Create a conda environment with the seq-seq-pan installation and activate it:

  ````sh
  conda create -n ssp seq-seq-pan

  # Activate the environment. To deactivate use "source deactivate"
  conda activate ssp 

  # update the snakemake version
  conda update -c conda-forge -c bioconda snakemake
  ````

### 4. Install Minimap2

  `````sh
  conda install minimap2
  ````
  
### 5. Install bedtools

  `````sh
  conda install bedtools
  ````

### 6. Install R and Rstudio

I you don'have it installed already:

R can be downloaded [here](https://cran.r-project.org/).
Rstudio can be downloaded [here](https://posit.co/download/rstudio-desktop/).

In R we will use the rtracklayer library to read the ATAC-seq bedgraph files. It would be possible to simply read in the ATAC-seq bedgraphs as a regular table, but because in reality you would probably want to work with compressed bigwig files, I show you how to use `rtracklayer`. This library is hosted by bioconductor and can be downloaded and installed in R as follows:

  ````r
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("rtracklayer")
  
  library(rtracklayer)
  ````

### 7. Download the data

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  To quickly run the entire code and have all the input files, clone [this](https://github.com/StevenVB12/Tutorial_pan_genomics) repository.
  
 ````sh
  git clone https://github.com/StevenVB12/Tutorial_pan_genomics
 ````
  
</div>
