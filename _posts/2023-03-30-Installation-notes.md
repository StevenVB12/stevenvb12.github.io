---
layout: post
title: "Installation notes"
categories: misc
---

### 1. Install Linux (Ubuntu) using WSL on windows (10 or higher)

(If you are on a Mac, just open a terminal and continue) 

* In your Windows search bar type `PowerShell`, right-click the result and "Run as administrator". 

* Enter the following command and restart your computer:
  ````sh
  wsl --install -d ubuntu
  ````

* Open Ubuntu by searching for it in your Windows search bar.

* Once Ubuntu has finished its initial setup you will need to create a username and password (this does not need to match your Windows user credentials).

* Finally, it’s always good practice to install the latest updates with the following commands, entering your password when prompted.
  ````sh
  sudo apt update
  sudo apt upgrade
  ````
  Press Y when prompted.
  
* To navigate to, for example, my google drive folder in Linux on windows, I do the following:
  ````sh
  cd /mnt/g/My Drive/Workshop_pan_genomics_12042023
  ````
  
### 2. Install Miniconda

> Conda is an open-source package management system and environment management system for installing and managing software packages and dependencies in various programming languages, including Python. It allows users to create and manage isolated environments with specific package versions to avoid conflicts between dependencies.
>
> Miniconda is a minimalistic version of the Anaconda (conda + Python/R) distribution. It includes only the essential components needed for creating a Python environment.

* First, download the installation [script](https://docs.conda.io/en/latest/miniconda.html#linux-installers).

* Run the installation script (change depending on your downloaded version):

  ````sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ````
  
* Type ‘yes’ to all questions.

* Close and reopen the shell.

### 3. Install seq-seq-pan 

#### Windows:

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
  Note that all subsequent programs will also be installed in the `ssp` environment. When you restart you terminal, you will need to reactivate it to run the programs.

#### Mac:

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #FFE3DD; border-color: #000000;">

  Sorry, due to a malfunctioning of the 'progressiveMauve' program the alignment step of seq-seq-pan doesn't work on Mac (unix), but you can dowload the outputs of the program [here](https://github.com/StevenVB12/Tutorial_pan_genomics/tree/main/seq-seq-pan_out).
  
</div>

To have the other functionalities of seq-seq-pan do the following:

* If you don't have git installed, check [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
* Clone the seq-seq-pan repository
  ````sh
  git clone https://gitlab.com/chrjan/seq-seq-pan.git
  cd seq-seq-pan
  chmod u+x seq-seq-pan*
  
  # Add seq-seq-pan to path so it is recognized as an executable from anywhere
  export PATH=/PathToFolder/seq-seq-pan/:$PATH
  ````

### 4. Install Minimap2

#### Windows:
  ````sh
  conda install minimap2
  ````
#### Mac:
  ````sh
  conda install "bioconda/label/cf201901" minimap2
  ````
  
### 5. Install bedtools

  ````sh
  conda install bedtools
  ````

### 6. Install R and Rstudio

I you don'have it installed already:

R can be downloaded [here](https://cran.r-project.org/).
Rstudio can be downloaded [here](https://posit.co/download/rstudio-desktop/).

> In R we will use the rtracklayer library to read the ATAC-seq bedgraph files. It would be possible to simply read in the ATAC-seq bedgraphs as a regular table, but because in reality you would probably want to work with compressed bigwig files, I show you how to use `rtracklayer`. This library is hosted by bioconductor and can be downloaded and installed in R as follows:

  ````r
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("rtracklayer")
  
  library(rtracklayer)
  ````

### 7. Download the data

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  To quickly run the entire code and have all the input files, clone [this](https://github.com/StevenVB12/Tutorial_pan_genomics) repository. If you don't have git installed, check [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
  
 ````sh
  git clone https://github.com/StevenVB12/Tutorial_pan_genomics
 ````
  
</div>
