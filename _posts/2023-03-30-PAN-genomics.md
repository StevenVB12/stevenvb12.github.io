---
layout: post
title: "Pan genome alignment tutorial"
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

In this tutorial we will align a piece of chromosome of two <i>Heliconius</i> butterfly species that includes the <i>optix</i> gene into a <strong>pan genome alignment</strong>. 

> <p align="center">
>  <img src="/docs/assets/Heliconius-melpomene-and-erato.jpg" width="300" title="rect()">
></p>
> Dorsal (top) and ventral (bottom) sides of <i>Heliconius melpomene rosina</i> (left) and <i>Heliconius erato demophoon</i> (right). The <i>optix</i> gene codes for a  transcription factor that plays a key role in the development of red color pattern elements. 

A pan genome refers to the complete set of genomic sequences shared by all individuals of a species, as well as the variable genomic sequences that are unique to specific individuals or subpopulations. We will use <strong>seq-seq-pan</strong> to construct the pan genome alignment, use some custom <strong>Python</strong> scripts to parse the outputs and use <strong>R</strong> to visualize the alignment. Additionally, we will transform the coordinates of <strong>Transposable Element (TE)</strong> annotations and <strong>chromatin accessibility</strong> profiles ([ATAC-seq](https://emea.illumina.com/techniques/popular-applications/epigenetics/atac-seq-chromatin-accessibility.html)) of developing head and wing tissues to the pan genome coordinate space and add them to this plot. 

This tutorial follows up on the figure we built in the [Minimap2 genome alignment tutorial](https://stevenvb12.github.io/misc/2023/03/30/Minimap2-alignment.html).

The final result should look like this:

<p align="center">
  <img src="/docs/assets/Plot_PAN.png" width="800" title="Minimap2">
</p>

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #3c763d; background-color: #dff0d8; border-color: #d6e9c6;">

  Note: I will sometimes switch between code used in the Linux terminal or Rstudio. I've put these in black and blue boxes, respectively.
  
</div>


### 2. Input data

(This is the same as in the [Minimap2 genome alignment tutorial](https://stevenvb12.github.io/misc/2023/03/30/Minimap2-alignment.html))

For this exercise, you can navigate to http://ensembl.lepbase.org/ and click on the links to the <i>Heliconius erato demophoon</i> (v1) and <i>Heliconius melpomene melpomene</i> (Hmel2) genome assemblies. At the right top, you can search for "optix". The search result will return you the gene model name (e.g. evm.TU.Herato1801.64) and its location (e.g. scaffold 'Herato1801' at position '1239943' to '1251211'). 
  
  ````
  # optix H. erato
  Herato1801:1239943-1251211

  # optix H. melpomene
  Hmel218003:705604-706407
  ````
  
When you click on the gene model, you can see what genes surround the <i>optix</i> gene. You can also see a botton on the left "Export data". This allows you to export the sequence as a `.fasta` file. Using this, you can try exporting not just the <i>optix</i> gene but a 2,000,000 bp region surrounding it.

You can also find these .fasta files [here](https://github.com/StevenVB12/Tutorial_pan_genomics/tree/main/sequences).
  ````
  # scaffold Herato1801 start 1 end 2000000
  sequences/Herato1801_1_2000000.fasta

  # scaffold Hmel213003 start 1 end 2000000
  sequences/Hmel213003_1_2000000.fasta
  ````


### 3. seq-seq-pan aligner

We can use [seq-seq-pan](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4401-3) to assemble the sequences in the fasta files into a pan genome. Seq-seq-pan extends the functionality of the multiple genome aligner [progressiveMauve](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892488/) by constructing a composite consensus or pan-genome that includes both homologous sequences or locally collinear blocks (LCBs) as well as lineage-specific (non-homologous) sequences in each of the genomes. This pan-genome is then used as the reference coordinates space for the multi genome alignment, which includes sequences specific to any of the genomes.

For our sequences, we will use seq-seq-pan as follows:

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````properties
  seq-seq-pan-wga --config genomefile=genome_list.txt outfilename=seq-seq-pan_out/SeqSeqPan_erato_melp_optix
  ````
  
</div>

> The `genome_list.txt` file includes a list (one per line) of the fasta sequences to be included in the pan genome assembly. The file can be downloaded [here](https://github.com/StevenVB12/Tutorial_pan_genomics/blob/main/genome_list.txt).

Seq-seq-pan will output several files. Two will be relevenat for us here:

> * The `_consensus.fasta` file includes a complete fasta sequence of the consensus pan genome (stitching all non-homologous sequences into the assembly and taking the allele that is most frequent among the multiple aligned genomes). This consensus file demarks the <strong>coordinate</strong> space of the pan genome and will be used when we want to map any positions in the original genomes (e.g. TE positions) to the pan genome.
> * The `.xmfa` file includes a list of the locally collinear blocks (LCBs). We will use this file to identify the sequences that are homologous or species-specific.
> * `> 1:` sequence identifier of the first genome in the `genome_list.txt` file.
> * `> 2:` sequence identifier of the second genome in the `genome_list.txt` file. (and so on)
> * `=` demarks separate LCBs.
> * `-` gaps in the aligned LCBs.

That's it, we have a pan genome!

### 4. Shared and unique sequences

We can now try to identify what parts of the sequences is identified as homologous or species-specific in the pan genome. We will use a custom Python script for this, available [here](https://github.com/StevenVB12/Tutorial_pan_genomics/blob/main/seq-seq-pan_blocks_intervals.py).

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````properties
  # For the Python script to work, we first need to remove newlines (enters) in the file from lines that include the sequence. This can be done with this perl oneliner:
  perl -pe 'chomp if /^[ATCGNSBDHVMRWYK-]/' seq-seq-pan_out/SeqSeqPan_erato_melp_optix.xmfa| sed 's/\=/\n\=/g' | sed 's/>/\n>/g' | sed '/^$/d' > seq-seq-pan_out/SeqSeqPan_erato_melp_optix.noNewline.xmfa

  # Now we can run the python script ('-g 1,2' is a list of the genome identifiers in the order of the genomes_list.txt file):
  python seq-seq-pan_blocks_intervals.py -I seq-seq-pan_out/SeqSeqPan_erato_melp_optix.noNewline.xmfa -g 1,2
  
  # The script annoyingly produces a 'TAB' at the end of each line which would break bedtools in the next step. We can remove that as follows:
  sed -i 's/[[:space:]]*$//' seq-seq-pan_out/1_blocks_intervals.bed 
  sed -i 's/[[:space:]]*$//' seq-seq-pan_out/2_blocks_intervals.bed 
  ````
  
</div>

> The output are files with `| chromosome | start | end |` positions of sequences in each genome but in the coordinate space of the pan genome (hence, a new line is generate when teh sequence is interrupted by a species-specific sequence in another genome). Such a format with start and end positions is typically called a `.bed` file format.

Next, we can use [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to subtract these files and get the unique sequences in each genome. 

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````properties
  bedtools subtract -sorted -a 1_blocks_intervals.bed -b 2_blocks_intervals.bed > blocks_unique_1.bed
  bedtools subtract -sorted -b 1_blocks_intervals.bed -a 2_blocks_intervals.bed > blocks_unique_2.bed
  ````
  
</div>

We also have a custom Python scripts to calculate sequence identity within the LCBs, available [here](https://github.com/StevenVB12/Tutorial_pan_genomics/blob/main/seq-seq-pan_bedfile_conservation.py).

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #000000; background-color: #000000; border-color: #000000;">
  
  ````properties
  # First, we need to transform the .xmfa file to a .fasta alignment
  python seq-seq-pan_toFasta.py -I seq-seq-pan_out/SeqSeqPan_erato_melp_optix.noNewline.xmfa -g 1,2
  cat genome* > seq-seq-pan_out/SeqSeqPan_erato_melp_optix.fasta
  
  # Next, we need to find the sequences that are shared between genomes. We can do this with bedtools intersect.
  bedtools intersect -sorted -a seq-seq-pan_out/1_blocks_intervals.bed -b seq-seq-pan_out/2_blocks_intervals.bed > seq-seq-pan_out/blocks_shared_1_2.bed
  
  # Finally, we can calculate the conservation. This will output a .bed like file with identity scores for each shared interval.
  python seq-seq-pan_bedfile_conservation.py -I seq-seq-pan_out/SeqSeqPan_erato_melp_optix.fasta -g 1,2 -b seq-seq-pan_out/blocks_shared_1_2.bed -o seq-seq-pan_out/conservation_shared_1_2.bed
  ````
  
</div>

### 5. Mapping annotations to the pan genome

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i optix_Hmel2.toMap.txt -n optix_Hmel2.toMap.pan

### 6. Visualizations in R

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i erato_5th_brain_start.toMap.txt -n erato_5th_brain_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i erato_5th_brain_end.toMap.txt -n erato_5th_brain_end.pan

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i erato_5th_FW_start.toMap.txt -n erato_5th_FW_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i erato_5th_FW_end.toMap.txt -n erato_5th_FW_end.pan

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i melp_5th_brain_start.toMap.txt -n melp_5th_brain_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i melp_5th_brain_end.toMap.txt -n melp_5th_brain_end.pan

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i melp_5th_FW_start.toMap.txt -n melp_5th_FW_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i melp_5th_FW_end.toMap.txt -n melp_5th_FW_end.pan


seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i miniMap_target_start.toMap.txt -n miniMap_target_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i miniMap_target_end.toMap.txt -n miniMap_target_end.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i miniMap_query_start.toMap.txt -n miniMap_query_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i miniMap_query_end.toMap.txt -n miniMap_query_end.pan

seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i TE_erato_start.toMap.txt -n TE_erato_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i TE_erato_end.toMap.txt -n TE_erato_end.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i TE_melp_start.toMap.txt -n TE_melp_start.pan
seq-seq-pan map -c SeqSeqPan_erato_melp_optix_consensus.fasta -p ./ -i TE_melp_end.toMap.txt -n TE_melp_end.pan


We're done! You should now see the figure that was at the top of this tutorial.
