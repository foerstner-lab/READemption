Performing an example analysis
==============================

Here you will be guided trough a small example analysis using a
publicly available RNA-Seq from NCBI GEO that was part of a
publication by `Kröger et
al. <http://www.ncbi.nlm.nih.gov/pubmed/24331466>`_. This is a
transcriptome analysis of *Salmonella* Typhimurium SL1344 in different
conditions. We will generate several output files in different
formats. The CSV (tabular separated plain text files) files can be
opened with any spreadsheet program like `LibreOffice
<https://www.libreoffice.org/>`_ or Excel. For inspecting the mappings
(in BAM format) and coverage files (wiggle format) you can use a
genome browser for example `IGB <http://bioviz.org/igb/>`_ or `IGV
<https://www.broadinstitute.org/igv/home>`_.

Generating a project
--------------------

At first we have to create the analysis folder and its subfolder. For
this we use the ``create`` subcommand::

  $ reademption create READemption_analysis
  Created folder "READemption_analysis2" and required subfolders.
  Please copy read files into folder "READemption_analysis2/input/reads" and reference sequences files into folder "READemption_analysis2/input/reference_sequences".

This will result in a folder structure as shown here:
::

 READemption_analysis
 ├── input
 │   ├── annotation_files
 │   ├── reads
 │   └── reference_sequences
 └── output
    ├── align
    │   ├── alignments
    │   ├── index
    │   ├── processed_reads
    │   ├── reports_and_stats
    │   │   ├── stats_data_json
    │   │   └── used_reademption_version.txt
    │   └── unaligned_reads
    ├── coverage
    │   ├── coverage-raw
    │   ├── coverage-tnoar_mil_normalized
    │   └── coverage-tnoar_min_normalized
    ├── deseq
    │   ├── deseq_raw
    │   └── deseq_with_annotations
    ├── gene_quanti
    │   ├── gene_quanti_combined
    │   └── gene_quanti_per_lib
    ├── viz_align
    ├── viz_deseq
    └── viz_gene_quanti


Retrieving the input data
-------------------------

We have to download the reference sequence (FASTA format) as well as
the annotation file (GFF3 format) for *Salmonella* from NCBI. As we
will use the URL of *Salmonella* Typhimurium SL1344's source FTP
folder it several times we store it in an environment variable called
``FTP_SOURCE``.  

::

  $ FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645

We download the reference sequence (the chromosome and three plasmids)
in FASTA format and store them in the ``reference_sequences``
folder. The files are saved with a different suffix (``.fa`` instead
of ``.fna``) as some genome browser (e.g. IGB) will not accept them as
FASTA files otherwise.

::
   
   $ wget -O READemption_analysis/input/reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
   $ wget -O READemption_analysis/input/reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
   $ wget -O READemption_analysis/input/reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
   $ wget -O READemption_analysis/input/reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna

We have to modify the header of the FASTA files as the sequence IDs
have to be the same as the ones in the first column of the GGF3 files
(see below) to be used in the gene quantification. This will be also
necessary if both, FASTA and GFF3 files, will be loaded in the IGB.

::

   $ sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/reference_sequences/NC_016810.fa
   $ sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/reference_sequences/NC_017718.fa
   $ sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/reference_sequences/NC_017719.fa
   $ sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/reference_sequences/NC_017720.fa

Then we download the GFF3 files that contain the annotations.
::

   $ wget -P READemption_analysis/input/annotations $FTP_SOURCE/*gff

Finally, we need the reads of the RNA-Seq libraries. To save some time
for running this examples we will work with subsampled libraries of 1M
reads each. This will the limit informative value of the results which
is acceptable as we just want to understand the workflow of the
READemption. Please be aware that READemption does not perform quality
trimming or adapter clipping so far. For this purpose use the `FASTX
toolkit <http://hannonlab.cshl.edu/fastx_toolkit/>`_, `cutadapt
<https://code.google.com/p/cutadapt/>`_ or other tools.

::

   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2

We have now all the necessary data available. The input folder should
look like this now:

::

   $ ls READemption_analysis/input/* 
   READemption_analysis/input/annotations:
   NC_016810.gff  NC_017718.gff  NC_017719.gff  NC_017720.gff
   
   READemption_analysis/input/reads:
   InSPI2_R1.fa.bz2  InSPI2_R2.fa.bz2  LSP_R1.fa.bz2  LSP_R2.fa.bz2
 
   READemption_analysis/input/reference_sequences:
   NC_016810.fa  NC_017718.fa  NC_017719.fa  NC_017720.fa

Processing and aligning the reads
---------------------------------

The first step it the read processing and mapping. Via parameters we
tell READemption to use 4 CPU (``-p 4``) and perform a poly-A-clipping
(``--poly_a_clipping``) before the mapping.

::

   $ reademption align -p 4 --poly_a_clipping READemption_analysis

Once this the mapping is done the file ``read_alignment_stats.csv`` is
created which can be found in
``READemption_analysis/output/align/reports_and_stats/``. It contains
several mapping statistics for example how many reads are successfully
aligned in total and how many were aligned to each replicon. We see
that more than 98 % of the reads are mapped for each library. Sorted
and indexed alignements in BAM format are stored in
``READemption_analysis/output/align/alignments``. We could load them
into a genome browser but instead we continue with the next step.


Generating coverage files
-------------------------

In order to generate strand specific coverage files with different
normalizations we use the subcommand ``coverage``.

::

   $ reademption coverage -p 4 READemption_analysis

The sets are stored in subfolder of
``READemption_analysis/output/coverage/``. The most oftenly used set
is stored in ``coverage-tnoar_min_normalized``. Here the coverage
values are normalized by the total number of aligned reads (TNOAR) of
the individual library and then multiplied by the lowest TNOAR value
of all libraries. These files could be inspected for differential
RNA-Seq (dRNA-Seq - comparing libraries with and without Terminator
Exonuclease treatment) data in order to determine transcriptional
start sites. They can be loaded in common genome browsers like `IGB
<http://bioviz.org/igb/>`_ or `IGV
<https://www.broadinstitute.org/igv/home>`_. Keep in mind that the
coverages of the reverse strand have negative values so you have to
adapt the scaling in some genome browsers.

Performing gene wise quantification
-----------------------------------

In this step we want to quantify the number of reads overlapping with
the locations of the annotation entries. With the ``--features``
parameter we configure ``reademption`` to just quantify CDS, tRNA and
rRNA entries.

::

   $ reademption gene_quanti -p 4 --features CDS,tRNA,rRNA READemption_analysis

After the quantification we find tables that contain the combined
counting for all entries in
``READemption_analysisoutput/gene_quanti/gene_quanti_combined/``. The
countings for mappings in sense and anti-sense are separately
listed. Besides the raw countings there are also tables for
countings normalized by the total number of reads and RPKM values.

Performing differential gene expression analysis
------------------------------------------------

To compare the gene expression of different conditions we apply the
subcommand ``deseq`` which makes use of the R library `DESeq2
<http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_. 

::

   $ reademption deseq \
      -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 \
      -c InSPI2,InSPI2,LSP,LSP READemption_analysis

::
  
We have to tell READemption which libraries are replicates of which
condition. This is done by the parameter ``-l`` and ``-c``. ``-l``
should hold a comma separated list of the libraries and ``-c`` the
corresponding conditions. In our case we have 4 libraries
(``InSPI2_R1.fa.bz2``, ``InSPI2_R2.fa.bz2``, ``LSP_R1.fa.bz2``,
``LSP_R2.fa.bz2``) and two condition (which we call ``InSPI2`` and
``LSP``). Just to make this association easier to understand:

::
   
      InSPI2_R1.fa.bz2  InSPI2_R2.fa.bz2  LSP_R1.fa.bz2  LSP_R2.fa.bz2 
         |                 |               |              |
      InSPI2            InSPI2            LSP            LSP 

When you call ``deseq`` it will compare all conditions with each other
and you can pick the comparison that you need. The raw ``DESeq2``
results are enriched with the original annotation information and are
stored in
``READemption_analysis/output/deseq/deseq_with_annotations/``

Create plots
------------

Finally we generate plots that visualize the results of the different
steps. ``viz_align`` creates histograms of the read length
distribution for the untreated and treated reads (saved in
``READemption_analysis/output/viz_align/``).

::
   
   $ reademption viz_align READemption_analysis

``viz_gene_quanti`` visualizes the gene wise countings. In our example
you will see that - as expected - the replicates are more similar to
each other than to the libs of the other condition. It also generates
bar plots that show the distribution of reads inside the different RNA
classes.

::

   $ reademption viz_gene_quanti READemption_analysis

``viz_deseq`` generates MA-plots as well as volcano plots.

::

   $ reademption viz_deseq READemption_analysis

