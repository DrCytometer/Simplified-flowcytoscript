# Simplified-flowcytoscript for Mac users
Mac-compatible version. For Windows, you'll be better off with the script on the Main branch.

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
* Mac-compatible. This version uses
* Speed. Optimizations throughout, although speed will vary depending on multithreading. This version will be slower than the Windows version, generally, because the requirement for OpenMP (difficult to install) was removed, so some parts are not parallelized. Phenograph in particular may be slow with large datasets.
* Both FCS and CSV files are accepted as input types.
* FCS data are automatically transformed as best befits the cytometer used. This avoids potentially serious issues with scaling of the data by inexperienced users. At present, only a few flow cytometers are supported: Aurora, ID7000, Fortessa, Symphony, Cytof (Helios) and FACSDiscoverS8. Note that S8 data is FCS3.2 format, which is not fully supported by FlowCore. I've not had issues apart from warning messages. If you are using a different cytometer, contact me. Transformation for Cytof data is Arcsinh while flow cytometry data will be with biexponential.
* Clustering can be performed using Phenograph or FlowSOM (via EmbedSOM).
* Clusters are automatically identified and named via matching to a cell type database. Take this with a grain of salt and check the results.
* tSNE performed in line with OptSNE modifications to learning rate.

Using the script: Install R and RStudio. Required packages will be installed/updated when the script is initiated.

In your favorite flow cytometry data analysis program (FlowJo, FCS Express), gate on the cells you wish to analyze and export those cells in new fcs or csv files. While exporting, adding group or variable tags to the file names will help you sort the files with the script.

To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". [More info](https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/)

Automated cluster identification adapted from [sc-type](https://github.com/IanevskiAleksandr/sc-type) by IanevskiAleksandr for scRNA-Seq data.

Create a folder for your analysis (preferably not in Dropbox or OneDrive). In this folder, put these items: A copy of flowcytoscript.Rmd A copy of 00_source_files Your files, inside a folder called “Data” Double click on the flowcytoscript file to open it in Rstudio. Run each code chunk in order by clicking on the green arrow in the upper right corner of the chunk.

Read through the presentation "Simplified flowcytoscript--instructions for use.pptx" for more detail, or watch the tutorial videos on [YouTube](https://www.youtube.com/watch?v=6x3Gwyf7-ww&t=3s).

For sample data and a demo analysis that you can try to recreate, visit [Dropbox](https://www.dropbox.com/scl/fo/s9h6z1k3rvliczv08uk6c/AGFFDoxnF1ttcZ7lTvddAQQ?rlkey=d3b224522jgq9g3rds8bnb3s9&dl=0).

Original publication on the [Crossentropy test](https://www.cell.com/cell-reports-methods/pdfExtended/S2667-2375(22)00295-8)

## Errors and bug reports
* Please save a copy of the notebook. This will produce an HTML document recording your entries and will facilitate troubleshooting.
