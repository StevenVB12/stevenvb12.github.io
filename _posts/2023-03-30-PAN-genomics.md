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






We're done! You should now see the figure that was at the top of this tutorial.
