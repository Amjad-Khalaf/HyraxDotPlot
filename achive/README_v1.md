# HyraxDotPlot
<p align="center">

<img src="https://github.com/Amjad-Khalaf/HyraxDotPlot/assets/92156267/2229e1d5-937e-4855-97bd-7c49b2385c73">

<p/>

## Table of Contents

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Base run](#1-base-run)
  - [2. Load feature bed files](#2-load-feature-bed-files)
  - [3. Load track bed files](#3-load-track-bed-files)
- [Acknowledgements](#acknowledgements)
- [Why hyraxes?](#why-hyraxes)

## About

HyraxDotPlot allows you to quickly and easily generate interactive dot plots from a [nucmer alignment](https://github.com/mummer4/mummer), with the options of loading bed files that outline the location of genomic features (presence/absence bed files, e.g. for genes of interest or telomere sequence motifs) as highlighted regions on the dot plot; and as tracks which run parallel to the axes their respective fastas are loaded on (quantitative bed files, e.g. for read depth along the genome). You can explore the full potential of HyraxDotPlot by opening the example output file `demo.html` in your browser.

<br>

<p align="center">

<img src="https://github.com/Amjad-Khalaf/HyraxDotPlot/assets/92156267/1ec63d9c-ce6a-4829-bb1e-5e7b86627e51">

<p/>


## Installation

To run HyraxDotPlot, you will need a virtual environment with `Python3` and the following python libraries installed:

```
holoviews
bokeh
matplotlib
```

Installation is simple, and running the tool requires the script `hyraxdotplot.py` only.

```
git clone https://github.com/Amjad-Khalaf/HyraxDotPlot.git
cd HyraxDotPlot
python hyraxdotplot.py -h ##check if it's running
```

## Usage

### 1. Base run

<p align="justify">

To run HyraxDotPlot with base features, generate a `nucmer coords file (with flags -T -l)` for your alignment and an index file (using `samtools faidx`) for each of your fasta files.

The base run will generate an interactive dot plot, with: 
1. the ability to zoom in down to single base pair resolution
2. the ability select sequences of interest by tapping them (using the "Tap" tool) and storing their ids in a downloabable list 
3. the ability to easily retrieve query sequence ids, subject sequence ids and match nucleotide identity % by hovering over any alignment of interest (using the "Hover" tool).

<p/>

```
python hyraxdotplot.py --coords_file $FILE --x_index_file $FILE --y_index_file $FILE

options:
  -h, --help            show this help message and exit
  --coords_file FILE    nucmer coordinates file generated with -T and -l flags; x axis fasta
                        should be the query sequence and y axis fasta should be the subject
                        sequence (i.e. nucmer x.fasta y.fasta)
  --x_index_file FILE   samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE   samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT     (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --plot_title STR      (optional) plot title; default is 'HyraxDotPlot'
  --output STR          (optional) output file name; default is hyraxdotplot.html
  --plot_width INT      (optional) plot width; default is 800
  --plot_height INT     (optional) plot height; default is 600
```

<br>

### 2. Load feature bed files

<p align="justify">
HyraxDotPlot also supports the loading of "qualitative" bed files that outline the location of genomic features (presence/absence bed files, e.g. for genes of interest or telomere sequence motifs), referred to as "feature bed files". These bed files will be loaded as highlighted regions on the dot plot, with a button to toggle them on or off, and a drop down menu to choose colours from. You may load a single bed file for one of the fastas, or a bed file for each fasta. Below is an example of a bed file which can be loaded in this manner.
<p/>

```
##the fourth column will be ignored for feature bed files

ptg000001l      980000  990000  309
ptg000005l      0       10000   270
ptg000009l      0       10000   202
```
<br>


```
python hyraxdotplot.py --coords_file $FILE --x_index_file $FILE --y_index_file $FILE
--x_feature_bed_file $FILE --y_feature_bed_file $FILE

options:
  -h, --help                show this help message and exit
  --coords_file FILE        nucmer coordinates file generated with -T and -l flags; x axis fasta
                            should be the query sequence and y axis fasta should be the subject
                            sequence (i.e. nucmer x.fasta y.fasta)
  --x_index_file FILE       samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE       samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT         (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --plot_title STR          (optional) plot title; default is 'HyraxDotPlot'
  --output STR              (optional) output file name; default is hyraxdotplot.html
  --plot_width INT          (optional) plot width; default is 800
  --plot_height INT         (optional) plot height; default is 600
  --x_feature_bed_file FILE (optional) qualitative bed file for x-axis fasta to be loaded as
                            highlighted regions on dot plot
  --y_feature_bed_file FILE (optional) qualitative bed file for y-axis fasta to be loaded as
                            highlighted regions on dot plot
```

### 3. Load track bed files

<p align="justify">

In addition to loading bed files as highlighted regions on the dot plot, HyraxDotPlot can also load "quantitative" bed files as tracks which run parallel to the axes their respective fastas are loaded on, referred to as "track bed files". You may load a single bed file for one of the fastas, or a bed file for each fasta. Below is an example of a bed file which can be loaded in this manner. You can load both feature and track bed files on the same plot. 
<p/>

```
##the fourth column is used for track bed files, and needs to be numerical

ptg000001l      980000  990000  309
ptg000005l      0       10000   270
ptg000009l      0       10000   202
```

<br>

```
python hyraxdotplot.py --coords_file $FILE --x_index_file $FILE --y_index_file $FILE
--x_track_bed_file $FILE --x_track_title $STR --x_track_feature_name $STR --x_track_colour $STR
--y_track_bed_file $FILE --y_track_title $STR --y_track_feature_name $STR --y_track_colour $STR

options:
  -h, --help                  show this help message and exit
  --coords_file FILE          nucmer coordinates file generated with -T and -l flags; x axis fasta
                              should be the query sequence and y axis fasta should be the subject
                              sequence (i.e. nucmer x.fasta y.fasta)
  --x_index_file FILE         samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE         samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT           (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --plot_title STR            (optional) plot title; default is 'HyraxDotPlot Interactive Dot Plot'
  --output STR                (optional) output file name; default is hyraxdotplot.html
  --plot_width INT            (optional) plot width; default is 800
  --plot_height INT           (optional) plot height; default is 600
  --x_track_bed_file FILE     (optional) quantitative bed file to be loaded as a track for
                              x-axis fasta, e.g. windowed coverage bed file
  --x_track_title STR         (optional) plot title for track loaded on x-axis fasta
  --x_track_feature_name STR  (optional) feature name for track loaded on x-axis fasta
  --x_track_colour STR        (optional) colour of track loaded on x-axis fasta; default is #56B4E9
  --y_track_bed_file FILE     (optional) quantitative bed file to be loaded as a track for
                              y-axis fasta, e.g. windowed coverage bed file
  --y_track_title STR         (optional) plot title for track loaded on y-axis fasta
  --y_track_feature_name STR  (optional) feature name for track loaded on y-axis fasta
  --y_track_colour STR        (optional) colour of track loaded on y-axis fasta; default is #f49ac2
```

## Acknowledgements

A huge thank you to Dr Claudia C Weber and Dr Ying Sims who have inspired many aspects of HyraxDotPlot and helped me debug it.

Cartoon image of hyrax was generated using https://www.bing.com/images/create and the font for the top part of the logo was taken from https://www.fontspace.com/category/cartoon.

## Why hyraxes?

You may have wondered at some point as you're using this tool or reading this page, where hyraxes are involved at all. They aren't. They're just cute and I like them. I hope you enjoy your plots!

<div align="center">
  
https://github.com/Amjad-Khalaf/HyraxDotPlot/assets/92156267/f20acd50-9cae-45eb-96ae-8394306f93f3

<div/>
