---
layout: post
title: "Minimap2 genome alignment tutorial"
categories: misc
---

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  To quickly run the entire code and have all the input files, clone [this](https://github.com/StevenVB12/Tutorial_pan_genomics) repository.
  
 ````
  git clone https://github.com/StevenVB12/Tutorial_pan_genomics
 ````
  
  To make sure you have everything installed and setup, check [here](https://stevenvb12.github.io/2023/03/30/Installation-notes.html).
  
</div>

### 1. Introduction

In this tutorial we will align a piece of chromosome of two <i>Heliconius</i> butterfly species that includes the <i>optix</i> gene. The <i>optix</i> gene codes for a  transcription factor that plays a key role in the development of red color pattern elements. We will use <strong>Minimap2</strong> to identify genomic intervals with significant resemblance between the two species and use <strong>R</strong> to visualize this alignment. Additionally, we will add <strong>Transposable Element (TE)</strong> annotations and <strong>chromatin accessibility</strong> profiles ([ATAC-seq](https://emea.illumina.com/techniques/popular-applications/epigenetics/atac-seq-chromatin-accessibility.html)) of developing head and wing tissues to this plot. 

The final result should look like this:

<p align="center">
  <img src="/docs/assets/Plot_minimap.png" width="800" title="Minimap2">
</p>

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  Note: I will often switch between code used in the Linux terminal or Rstudio. I've put these in black and blue boxes, respectively.
  
</div>


### 2. Input data
  
For this exercise, you can navigate to http://ensembl.lepbase.org/ and click on the links to the <i>Heliconius erato demophoon</i> (v1) and <i>Heliconius melpomene melpomene</i> (Hmel2) genome assemblies. At the right top, you can search for "optix". The search result will return you the gene model name (e.g. evm.TU.Herato1801.64) and its location (e.g. scaffold 'Herato1801' at position '1239943' to '1251211'). 
  
  ````
  # optix H. erato
  Herato1801:1239943-1251211

  # optix H. melpomene
  Hmel218003:705604-706407
  ````
  
When you click on the gene model, you can see what genes surround the <i>optix</i> gene. You can also see a botton on the left "Export data". This allows you to export the sequence as a `.fasta` file. Using this, you can try exporting not just the <i>optix</i> gene but a 2,000,000 bp region surrounding it.

You can also find these .fasta files [here](https://github.com/StevenVB12/Tutorial_pan_genomics/tree/main/input).
  ````
  # scaffold Herato1801 start 1 end 2000000
  Herato1801_1_2000000.fasta.gz

  # scaffold Hmel213003 start 1 end 2000000
  Hmel213003_1_2000000.fasta.gz
  ````
  
  
### 3. Minimap2 aligner

We can use [Minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) to align the sequences in the fasta files. Minimap2 is intended to efficiently map error prone long-read sequences to a reference (e.g. from nanopore or Pacbio sequencing) but it also works well for pairwise alignment of whole chromosomes and genomes. 

For our sequences, we will use Minimap2 as follows:

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````
  minimap2 -x asm20 Hmel213003_1_2000000.fasta Herato1801_1_2000000.fasta --no-long-join -r 200 | cut -f 1-12 > Minimap_erato_melp.sam
  ````
</div>

Note that we have set here `Hmel213003_1_2000000.fasta` as the <strong>target</strong> and `Herato1801_1_2000000.fasta` as the <strong>query</strong> sequence.

The `-x asm20` option is a preset to optimize alignment of sequences divergent up to 20%. The `--no-long-join -r 200` option assures that alignments include no gaps longer than 200 bp (otherwise the result may be one long contiguous alignment). See [here](https://lh3.github.io/minimap2/minimap2.html) for more detailed information on the Minimap2 command. 

The `| cut -f 1-12` will only take the first twelve columns of the Minimap2 output. This is sent to the file `Minimap_erato_melp.sam`.

The output should look like this (without the column titles as shown here. See [here](https://lh3.github.io/minimap2/minimap2.html) for more info):

| queryName | queryLength | queryStart | queryEnd | char | targetName | targetLength | targetStart | targetEnd | matchingBases | matchLength | matchQuality |
| :------------- | :-------- | ------: | ------: | ---: | :------------- | :-------- | ------: | ------: | ------: | ------: | ------: |
| Herato1801 | 2000000 | 700514 | 710441 | + | Hmel218003 | 1895173 | 228710 | 238849 | 3366 | 10410 | 60 |
| ... | ... | ...  | ... | ... | ... | ... | ... | ... | ... | ... | ... |
  
Note here that the target length isn't actually 2000000 bp but 1895173 bp. That's because the <i>H. melpomene</i> scaffold Hmel218003 with the <i>optix</i> gene has a length of 1895173 bp only. 

That's it. For me the alignment step took less than a second. Thanks [Heng Li](https://en.wikipedia.org/wiki/Heng_Li)!
  
  
### 4. Visualizations in R

Now we can switch to Rstudio (but we will be generating some extra input files in Linux throughout the pipeline).

#### 4.1. Set you working directory:

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #31708f; background-color: #d9edf7; border-color: #bce8f1;">
  
  ```r
  setwd("G:/My Drive/Workshop_pan_genomics_12042023")
  ```

</div>

(You can also use the '...' under the 'File' tab on the right of Rstudio to navigate to your folder then click 'More' > 'Set as working directory'.

#### 4.2. Load the Minimap2 output:

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #31708f; background-color: #d9edf7; border-color: #bce8f1;">
  
  ```r
  miniMap_out <- read.table('Minimap_melp_erato.sam', header = FALSE, sep = '\t')
                 
  # set the column names
  colnames(miniMap_out) <- c('queryName', 'queryLength', 'queryStart', 'queryEnd', 'char', 'targetName', 'targetLength', 'targetStart', 'targetEnd', 'matchingBases', 'matchLength', 'matchQuality')
  
  # check the table
  head(miniMap_out)
  ```

</div>

#### 4.3. Plot the first minimap2 match

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #31708f; background-color: #d9edf7; border-color: #bce8f1;">
  
  ```r
  # First define x-axis coordinates. This will help us later define the genomic window (xlim) we want to zoom in on.
  start = 0
  end = 3000000

  # Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
  # For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
  plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

  # Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
  rect(0,8,2000000,9, col = 'deepskyblue', border = NA)
  rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

  # We also know the positions of the optix gene and we can add these with the same rect() function trick.
  rect(705604,9,706407,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

  rect(1239943,0,1239972,1, col = 'red', border = 'red') # first exon position optix erato
  rect(1250591,0,1251211,1, col = 'red', border = 'red') # second exon position optix erato
  rect(1239943,0.5,1251211,0.5, col = 'red', border = 'red') # A little line between the two and we have a gene model!

  # With the text function, we can add the gene and species names at the appropriate coordinates.
  text(706407, 9.5, substitute(paste(italic('optix'))), pos = 4)
  text(1251211, 9.5, substitute(paste(italic('optix'))), pos = 4)

  text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
  text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)
  ```

</div>
