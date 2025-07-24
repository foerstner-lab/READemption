Visualization of RNA-seq data
==============================

This guide explains how to visualize RNA-seq data files in integrated genome browser (IGB). 

Installing IGB
--------------

1. Install IGB for your system from IGB official website: https://www.bioviz.org/ .
2. Follow the installation requirements for Linux, Windows and macOS.
3. Launch IGB.

For linux, run the following command in terminal after installation to open IGB:
::

    $ chmod +x IGB-linux-amd64-<version>.sh
    $ ./IGB-linux-amd64-<version>.sh


Preparing data for visualization
--------------------------------

After running the READemption analysis, you"ll get wiggles files in coverage folder. 
These wiggle files will be used for visualization. 

1. Go to the coverage folder first.
2. Navigate ".wig" or ".wig.gz" files of your sample.
3. Unzip the .wig.gz files (optional) by running the following command:
         $ Gunzip *.wig.gz

Load the data in IGB
--------------------

1. Open IGB.
2. Load reference genome file (".fa", ".fna", ".fasta") and annotation file (".gff3") by clicking on
 "File" → "Open file" on the top left corner, and then click on "Load sequence" on top right corner next to zoom option .
3. Load the coverage ".wig" files. 
4. Zoom in and zoom out into the regions of your interest.

NOTE
====

After running `READemption_analysis`, you will get two folders of coverage:

- `READemption_analysis/output/salmonella_coverage-tnoar_mil_normalized/`
- `READemption_analysis/output/salmonella_coverage-tnoar_min_normalized/`

You can use either of these folders for visualization, but for IGB display, the following is **recommended** for a clearer view of normalized expressions:

- `READemption_analysis/output/salmonella_coverage-tnoar_mil_normalized/`



Tips for Visualization in IGB
=============================

1. You can change the color of the strands of the reference genome, GFF3, or wiggle file data by locating the track in the IGB main window where there is a table-like view named **Data Management View**. Click on the color appearing in the box below **FG**.

2. IGB only supports file formats like `.fasta`, `.gff3`, `.bigwig`, `.bam`, `.bed`, `.wig`, and `.vcf`. Make sure you load the correct file format.

3. Verify your files are compatible with your version of IGB.


