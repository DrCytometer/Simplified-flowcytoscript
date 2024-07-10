# Simplified-flowcytoscript
This version works best for Windows. For Mac, see the MacVersion branch.

A simplified, semi-automated high parameter workflow for people who don't know how to code. High dimensional flow cytometry analysis.
A simplified complete workflow in R for analysis of high parameter flow cytometry data including the Crossentropy method.

Code is available and free for academic users. Commercial users should contact Adrian Liston to discuss licensing options.

This simplified version of the flowcytoscript (Crossentropy test) is intended to be usable by people with little to no experience in programming. All inputs are via plain language prompts in an RStudio markdown notebook. All outputs are organized in folders as before, but additionally an HTML document summarizing the results is created with each analysis. Check out "flowcytoscript.nb.html" for an example of the output.

## Features
* Clustering with histograms of expression, barcharts and tables of frequencies
* Automated cluster identification
* tSNE and UMAP
* PCA on MFIs to show sample-level differences
* Crossentropy testing on tSNE and UMAP
* Statistical testing on markers and clusters
* Heatmaps, dendrograms

## Improvements
* Speed. Optimizations throughout should render this approximately 10x faster, although this will vary depending on multithreading.
* Both FCS and CSV files are accepted as input types.
* FCS data are automatically transformed as best befits the cytometer used. This avoids potentially serious issues with scaling of the data by inexperienced users. At present, only a few flow cytometers are supported: Aurora, ID7000, FACSDiscoverS8, Bio-Rad ZE5 (Yeti), Fortessa, Symphony. If you are using a different cytometer, contact me.
* Clustering can be performed using Phenograph or FlowSOM (via EmbedSOM).
* Clusters are automatically identified and named via matching to a cell type database. Take this with a grain of salt and check the results.
* tSNE performed in line with OptSNE modifications to learning rate.

Using the script: Install R, Rstudio and Rtools. To use this version with Mac, you’ll need command line tools and OpenMP. The flowcytoscript_setup.r script (in 00_source_files) can be used to facilitate set-up for new users of R. Installation for Mac is difficult due to OpenMP and the Mac-specific version (see branch MacVersion) avoids this, with some loss in speed.

In your favorite flow cytometry data analysis program (FlowJo, FCS Express), gate on the cells you wish to analyze and export those cells in new fcs or csv files. While exporting, adding group or variable tags to the file names will help you sort the files with the script.

To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". [More info](https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/)

Automated cluster identification adapted from [sc-type](https://github.com/IanevskiAleksandr/sc-type) by IanevskiAleksandr for scRNA-Seq data.

Create a folder for your analysis (preferably not in Dropbox or OneDrive). In this folder, put these items: A copy of flowcytoscript.Rmd A copy of 00_source_files Your files, inside a folder called “Data” Double click on the flowcytoscript file to open it in Rstudio. Run each code chunk in order by clicking on the green arrow in the upper right corner of the chunk.

Read through the presentation "Simplified flowcytoscript--instructions for use.pptx" for more detail, or watch the tutorial videos on [YouTube](https://www.youtube.com/watch?v=6x3Gwyf7-ww&t=3s).

For sample data and a demo analysis that you can try to recreate, visit [Dropbox](https://www.dropbox.com/scl/fo/s9h6z1k3rvliczv08uk6c/AGFFDoxnF1ttcZ7lTvddAQQ?rlkey=d3b224522jgq9g3rds8bnb3s9&dl=0).

Original publication on the [Crossentropy test](https://www.cell.com/cell-reports-methods/pdfExtended/S2667-2375(22)00295-8)

## Transformation details
For CyTOF (Helios) data, ArcSinh with a co-factor of 5 is used. 
```arcsinhTransform(a=1, b=1/5, c=0)```

For flow cytometry data, a FlowJo biexponential is used with the following co-efficients:

Aurora:
```
width.basis <- -1000
max.value <- 4194303
log.decades <- 5.5
```
ID7000:
```
width.basis <- -500
max.value <- 1000000
log.decades <- 5
```
FACSDiscoverS8:
```
width.basis <- -1000
max.value <- 3162277.6602
log.decades <- 5
```
ZE5:
```
width.basis <- -50
max.value <- 262144
log.decades <- 4.42
```
Others are transformed with coefficients appropriate for the BD Fortessa or Symphony:
```
width.basis <- -100
max.value <- 262144
log.decades <- 4.5
```

These values are selected as they represent more or less a default view of the data on the cytometer. The appropriateness
of these transforms assumes you've optimized your panel by visually inspecting the staining and adjusting marker
intensity/separation based on what you can see. These transforms will be best for panels that are relatively large for
the cytometer.

## FACSDiscoverS8 notes
The BD FACSDiscoverS8 produces FCS3.2 files, and FlowCore does not yet provide full support for this format. You will
get warning messages to this effect. I've yet to experience any actual problems.

The S8 files contain hundreds of parameters, most of which are probably irrelevant to your analysis. To simplify things,
the script automatically selects the SpectralFX-unmixed channels, discarding the raw, OLS unmixed and imaging parameters. 
As such, you can export all channels or all compensated channels from FlowJo. Some of the marker names (if they contain spaces)
will be truncated, though, and you'll have to edit these. If you have multiple channels that truncate to become identical (e.g.,
TCR beta -> TCR and TCR gd -> TCR), these will be listed in the same order they appear in the FCS file in your analysis software.

## Errors and bug reports
* Please save a copy of the notebook. This will produce an HTML document recording your entries and will facilitate troubleshooting.
* If you've run the script and encountered errors, you will likely need to start the analysis over in a new folder to
* avoid error propagation. Alternatively, delete the cache files (tsne_cache, umap_cache).
