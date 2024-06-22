# Simplified-flowcytoscript
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
* FCS data are automatically transformed as best befits the cytometer used. This avoids potentially serious issues with scaling of the data by inexperienced users. At present, only a few flow cytometers are supported: Aurora, ID7000, Fortessa, Symphony. If you are using a different cytometer, contact me.
* Clustering can be performed using Phenograph or FlowSOM (via EmbedSOM).
* Clusters are automatically identified and named via matching to a cell type database. Take this with a grain of salt and check the results.
* tSNE performed in line with OptSNE modifications to learning rate.

Using the script: Install R, Rstudio and Rtools. For Mac, you’ll need command line tools and OpenMP. The flowcytoscript_setup.r script (in 00_source_files) can be used to facilitate set-up for new users of R. Installation for Mac is difficult due to OpenMP and a Mac-specific version is in progress.

In your favorite flow cytometry data analysis program (FlowJo, FCS Express), gate on the cells you wish to analyze and export those cells in new fcs or csv files. While exporting, adding group or variable tags to the file names will help you sort the files with the script.

To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". [More info](https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/)

Automated cluster identification adapted from [sc-type](https://github.com/IanevskiAleksandr/sc-type) by IanevskiAleksandr for scRNA-Seq data.

Create a folder for your analysis (preferably not in Dropbox or OneDrive). In this folder, put these items: A copy of flowcytoscript.Rmd A copy of 00_source_files Your files, inside a folder called “Data” Double click on the flowcytoscript file to open it in Rstudio. Run each code chunk in order by clicking on the green arrow in the upper right corner of the chunk.

Read through the presentation "Simplified flowcytoscript--instructions for use.pptx" for more detail.

Original publication on the [Crossentropy test](https://www.cell.com/cell-reports-methods/pdfExtended/S2667-2375(22)00295-8)

## Errors and bug reports
* Please save a copy of the notebook. This will produce an HTML document recording your entries and will facilitate troubleshooting.
