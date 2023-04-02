---
layout: post
title: "Minimap2 genome alignment tutorial"
categories: misc
---

### 1. Introduction

In this tutorial we will align a piece of chromosome of two <i>Heliconius</i> butterfly species that includes the <i>optix</i> gene. The <i>optix</i> gene codes for a   transcription factor that plays a key role in the development of red color pattern elements. We will use <strong>Minimap2</strong> to identify genomic intervals with significant resemblance between the two species and use <strong>R</strong> to visualize this alignment. Additionally, we will add <strong>Transposable Element (TE)</strong> annotations and <strong>chromatin accessibility</strong> profiles ([ATAC-seq](https://emea.illumina.com/techniques/popular-applications/epigenetics/atac-seq-chromatin-accessibility.html)) of developing head and wing tissues to this plot. 

The final result should look like this:

<p align="center">
  <img src="/docs/assets/Plot_minimap.png" width="800" title="Minimap2">
</p>

### 2. Input data
  
  For this exercise, you can navigate to http://ensembl.lepbase.org/ and click on the links to the <i>Heliconius erato demophoon</i> (v1) and <i>Heliconius melpomene melpomene</i> (Hmel2) genomes. At the right top, you can search for "optix". The search result will return you the gene model name (e.g. evm.TU.Herato1801.64) and its location (e.g. scaffold 'Herato1801' at position '1239943' to '1251211'). 
  
  ````
    # optix H. erato
    Herato1801:1239943-1251211
    
    # optix H. melpomene
    Hmel218003:705604-706407
  ````
  
  When you click on the gene model, you can see what genes surround the <i>optix</i>. You can also see a botton on the left "Export data". This allows you to export the sequence as a .fasta file. Using this, you can try exporting not just the <i>optix</i> but a 2,000,000 bp region surrounding it.
  
  ````
  # scaffold Herato1801 start 1 end 2000000
  Herato1801_1_2000000.fasta.gz
  
  # scaffold Hmel213003 start 1 end 2000000
  Hmel213003_1_2000000.fasta.gz
  ````
  
  You can also find these fasta files [here](https://github.com/StevenVB12/Tutorial_pan_genomics/tree/main/input).
  
  ### 3. Minimap2 aligner
  
  The `-x asm20` option is a preset to optimize alignment of sequences divergent up to 20%. The `--no-long-join -r 200` option assures that long gaps aren't ignored for building long alignments and that alignments include no gaps longer than 200 bp. See [here](https://lh3.github.io/minimap2/minimap2.html) for more detailed information on the Minimap2 command. 
  
<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````
  minimap2 -x asm20 Hmel213003_1_2000000.fasta Herato1801_1_2000000.fasta --no-long-join -r 200 | cut -f 1-12 > Minimap_erato_melp.sam
  ````
  
</div>
  
  Note that we have set here `Hmel213003_1_2000000.fasta` as the <strong>target</strong> and `Herato1801_1_2000000.fasta` as the <strong>query</strong> sequence.
  
  ### 4. Visualizations in R
  
  ```r
  # First set you working directory.
  # (You can also use the '...' under the 'File' tab on the right of Rstudio 
  # to navigate to your folder then click 'More' > 'Set as working directory')

  setwd("G:/My Drive/Workshop_pan_genomics_12042023")
  ```















<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #31708f; background-color: #d9edf7; border-color: #bce8f1;">

  ````
  Note
  ````

  
</div>


<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  ````
  Note
  ````
  
</div>

