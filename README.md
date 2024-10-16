# HyraxDotPlot

<p align="center">

<img src="https://github.com/user-attachments/assets/8219b9aa-c78c-4c68-a082-96590c80876e" width="800">

</p>

## About

<div align="justify">

**H**<sub>(tml-based-output)</sub>yrax**DotPlot** is a tool written in Python (using the library [bokeh](https://docs.bokeh.org/en/latest/), which allows you to quickly and easily generate interactive dot plots from a pair-wise genome file (a [minimap2](https://github.com/lh3/minimap2) paf file, a [FASTGA](https://github.com/thegenemyers/FASTGA) paf file, or a [nucmer coords alignment](https://github.com/mummer4/mummer). It provides you with the options below.

1. Loading **quantitative** genomic feature bed files (e.g. for read depth along the genome) as tracks which run parallel to the axes their respective fastas are loaded on.
2. Loading **qualitative** genomic feature bed files (presence/absence bed files, e.g. for genes of interest or telomere sequence motifs) as highlighted regions on the dot plot.
3. Loading **qualitative** genomic feature bed files (e.g. for gene annotation bed files) as tracks which run parallel to the axes their respective fastas are loaded on.

You can load feature bed files, track bed files, and annotation track bed files on the same plot (or in any combination desired).

HyraxDotPlot is lightweight, and produces a static png, a static svg, and a interactive html file which can be opened in a browser of your choice. In the html output file, you can zoom in down to 1 base-pair resolution, and if the flag `--curation_mode` is used, you can select sequences of interest by tapping them (using the "Tap" tool) and storing their ids in a downloabable list. Additionally, using the "Hover" tool, you also have the ability to easily retrieve query sequence ids, subject sequence ids and match nucleotide identity % by hovering over any alignment of interest. 

## Update: What's new in version 2?

1. paf files can now be used, as long as they contain the `de` tag.
2. The `--size_threshold` parameter allows you to control the minimum alignment size to be plotted. This is crucial if your alignment of choice does not carry out any chaining, or you are plotting the alignment of two large genomes. I recommend using a minimum of `--size_threshold 1000` for rapid execution (this is set by default).
3. `--x_annotation_bed_file` and `--y_annotation_bed_file` allow you to load gene/repeat/other annotation bed files as tracks which run parallel to the axes their respective fastas are loaded on.
4. A colour bar is added at the bottom of all plots automatically to show the identity of the different alignments plotted.
5. The ability to select contigs/scaffolds/chromosomes/sequences and add them to a downloadable list used to be automatic in the previous version, but now can only be enabled using `--curation_mode`. This was done to increase efficiency.
6. By default, both a static png and a static svg are generated in addition to the interactive html output file. Even if `--curation_mode` is enabled, the relevant selection and download widgets are not displayed in the static outputs.

</div>

## Installation

<div align="justify">
  
To run HyraxDotPlot, you will need a virtual environment with `Python3` and the following installed.

</div>

```
holoviews==1.19.1 
bokeh==3.6.0
matplotlib==3.9.2
selenium==4.25.0
chromedriver==2.24.1
webdriver-manager==4.0.2
firefox==131.0 
geckodriver==0.35.0 
```

<div align="justify">
  
Installation is simple, and running the tool requires the script `hyraxdotplot.py` only.

</div>

```
git clone https://github.com/Amjad-Khalaf/HyraxDotPlot.git
cd HyraxDotPlot
python hyraxdotplot.py -h ##check if it's running
```

## Usage

### 1. Base run

<div align="justify">

To run HyraxDotPlot with base features, generate a `nucmer coords file (with flags -T -l)` or a `paf file with the de tag` for your alignment and an index file (using `samtools faidx`) for each of your fasta files.

</div>

```
python hyraxdotplot.py --coords_file $FILE --x_index_file $FILE --y_index_file $FILE

options:
  -h, --help             show this help message and exit
  --coords_file FILE     nucmer coordinates file generated with -T and -l flags; x axis fasta
                         should be the query sequence and y axis fasta should be the subject
                         sequence (i.e. nucmer x.fasta y.fasta)
  --paf_file FILE        minimap2 paf file generated with -c flag or FastGA paf file with the
                         "de" tag available; x-axis fasta should be the query sequence and y-axis
                         fasta should be the subject sequence (i.e. minimap2/FastGA x.fasta y.fasta)
  --x_index_file FILE    samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE    samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT      (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --size_threshold FLOAT (optional) minimum length of matches to be plotted; default is 1000
  --plot_title STR       (optional) plot title; default is 'HyraxDotPlot'
  --output_prefix STR    (optional) output file name prefix; default is hyraxdotplot
  --plot_width INT       (optional) plot width; default is 800
  --plot_height INT      (optional) plot height; default is 600
  --curation_mode        (optional) flag that adds tap tool and allows you to select sequences of interest
```

<p align="center">
<img src="https://github.com/user-attachments/assets/8f817163-162d-4d5c-ad0f-f061a5c37a13">
</p>


### 2. Load feature bed files

<div align="justify">
  
HyraxDotPlot also supports the loading of "qualitative" bed files that outline the location of genomic features (presence/absence bed files, e.g. for genes of interest or telomere sequence motifs), referred to as "feature bed files". These bed files will be loaded as highlighted regions on the dot plot, with a button to toggle them on or off, and a drop down menu to choose colours from. You may load a single bed file for one of the fastas, or a bed file for each fasta. Below is an example of a bed file which can be loaded in this manner. You can load feature bed files, track bed files, and annotation track bed files on the same plot (or in any combination desired).

</div>

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
  --paf_file FILE           minimap2 paf file generated with -c flag or FastGA paf file with the
                            "de" tag available; x-axis fasta should be the query sequence and y-axis
                            fasta should be the subject sequence (i.e. minimap2/FastGA x.fasta y.fasta)
  --x_index_file FILE       samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE       samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT         (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --size_threshold FLOAT    (optional) minimum length of matches to be plotted; default is 1000
  --plot_title STR          (optional) plot title; default is 'HyraxDotPlot'
  --output_prefix STR       (optional) output file name prefix; default is hyraxdotplot
  --plot_width INT          (optional) plot width; default is 800
  --plot_height INT         (optional) plot height; default is 600
  --x_feature_bed_file FILE (optional) qualitative bed file for x-axis fasta to be loaded as
                            highlighted regions on dot plot
  --y_feature_bed_file FILE (optional) qualitative bed file for y-axis fasta to be loaded as
                            highlighted regions on dot plot
  --curation_mode           (optional) flag that adds tap tool and allows you to select sequences of interest
```

<p align="center">
<img src="https://github.com/user-attachments/assets/4639984e-d080-4a25-ba06-2b6e945d3324">
</p>

### 3. Load track bed files

<div align="justify">

In addition to loading bed files as highlighted regions on the dot plot, HyraxDotPlot can also load "quantitative" bed files as tracks which run parallel to the axes their respective fastas are loaded on, referred to as "track bed files". You may load a single bed file for one of the fastas, or a bed file for each fasta. Below is an example of a bed file which can be loaded in this manner. You can load feature bed files, track bed files, and annotation track bed files on the same plot (or in any combination desired). 

</div>

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
  --paf_file FILE             minimap2 paf file generated with -c flag or FastGA paf file with the
                              "de" tag available; x-axis fasta should be the query sequence and y-axis
                              fasta should be the subject sequence (i.e. minimap2/FastGA x.fasta y.fasta)
  --x_index_file FILE         samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE         samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT           (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --size_threshold FLOAT      (optional) minimum length of matches to be plotted; default is 1000
  --plot_title STR            (optional) plot title; default is 'HyraxDotPlot Interactive Dot Plot'
  --output_prefix STR         (optional) output file name prefix; default is hyraxdotplot
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
  --curation_mode             (optional) flag that adds tap tool and allows you to select sequences of interest
```

<p align="center">
<img src="https://github.com/user-attachments/assets/30d4fea5-d999-46a1-b319-a865e5a70e19">
</p>

### 4. Load annotation track bed files

<div align="justify">
  
HyraxDotPlot can also load "qualitative" bed files as tracks which run parallel to the axes their respective fastas are loaded on, referred to as "annotation bed files". You may load a single bed file for one of the fastas, or a bed file for each fasta. Below is an example of a bed file which can be loaded in this manner. You can load feature bed files, track bed files, and annotation track bed files on the same plot (or in any combination desired).

</div>

```
##you need: column1=sequence_id, column2=start, column3=end, column4=feature_id, column6=strand

ptg000024l      110474  120107  1at6029 1171.7  +
ptg000055l      102871  110374  11at6029        694.6   +
ptg000634l      20162   25100   14at6029        730.3   +
```
<br>

```
python hyraxdotplot.py --coords_file $FILE --x_index_file $FILE --y_index_file $FILE
--x_annotation_bed_file $FILE --y_annotation_bed_file $FILE

options:
  -h, --help                    show this help message and exit
  --coords_file FILE            nucmer coordinates file generated with -T and -l flags; x axis fasta
                                should be the query sequence and y axis fasta should be the subject
                                sequence (i.e. nucmer x.fasta y.fasta)
  --paf_file FILE               minimap2 paf file generated with -c flag or FastGA paf file with the
                                "de" tag available; x-axis fasta should be the query sequence and y-axis
                                fasta should be the subject sequence (i.e. minimap2/FastGA x.fasta y.fasta)
  --x_index_file FILE           samtools .fai file of fasta to be on x-axis of dotplot
  --y_index_file FILE           samtools .fai file of fasta to be on y-axis of dotplot
  --threshold FLOAT             (optional) minimum nucleotide identity of matches to be plotted; default is 90
  --size_threshold FLOAT        (optional) minimum length of matches to be plotted; default is 1000
  --plot_title STR              (optional) plot title; default is 'HyraxDotPlot Interactive Dot Plot'
  --output_prefix STR           (optional) output file name prefix; default is hyraxdotplot
  --plot_width INT              (optional) plot width; default is 800
  --plot_height INT             (optional) plot height; default is 600
  --x_annotation_bed_file FILE  (optional) gene or repeat annotation bed file for x-axis fasta
  --y_annotation_bed_file FILE  (optional) gene or repeat annotation bed file for y-axis fasta
  --curation_mode               (optional) flag that adds tap tool and allows you to select sequences of interest
```

<p align="center">
<img src="https://github.com/user-attachments/assets/ce525138-92e7-4216-a393-fa3d19a5c17c">
</p>

## Acknowledgements

<div align="justify">
  
A huge thank you to Dr Claudia C Weber and Dr Ying Sims who have inspired many aspects of HyraxDotPlot and helped me debug it.

</div>

## Why hyraxes?

<div align="justify">
  
You may have wondered at some point as you're using this tool or reading this page, where hyraxes are involved at all. They aren't. They're just cute and I like them. I hope you enjoy your plots!

</div>

<div align="center">
  
https://github.com/Amjad-Khalaf/HyraxDotPlot/assets/92156267/f20acd50-9cae-45eb-96ae-8394306f93f3

</div>
