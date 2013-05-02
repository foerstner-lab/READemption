Performing a standard analysis
==============================

In the following we will run a small example analysis using a publicly
available RNA-Seq data set. 

Generating a project
--------------------

Creating a new project::

  $ python3.3 rapl/rapl.py create my_rna_seq_analysis
  Created folder "my_rna_seq_analysis" and required subfolders.
  Please copy read files into folder "my_rna_seq_analysis/input/reads" and reference sequences files into folder "my_rna_seq_analysis/input/reference_sequences".


  $ ls my_rna_seq_analysis/*
  my_rna_seq_analysis/input:
  annotation_files  reads  reference_sequences

  my_rna_seq_analysis/output:
  coverages-raw                   deseq_comparisons           read_alignments-index            reports_and_stats
  coverages-tnoar_mil_normalized  gene_wise_quantifications   read_alignments-processed_reads  stats_data_json
  coverages-tnoar_min_normalized  read_alignments-alignments  read_alignments-unaligned_reads

Retrieving the input data
-------------------------


Store the URL of source FTP folder in an environment variable.
::
   $ FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645

Download the reference sequence (the chromosome and three plasmids) in Fasta format.
::
   $ wget -cP my_rna_seq_analysis/input/reference_sequences ${FTP_SOURCE}/NC_016810.fna ...

Download the annotation for in GFF format.
::
   $ wget -cP ${FTP_SOURCE}/NC_016810.gff ...

Download the reads in FASTA format
::
   $ wget ${FTP_SOURCE}/NC_016810.fna ...

Aligning the reads to the reference genome
------------------------------------------

Generating coverage files
-------------------------

`Integrated genome browser (IGB) <http://bioviz.org/>`_ 
`Integrative genome viewer (IGV) <https://www.broadinstitute.org/software/igv/>`_


Performing gene wise quantification
-----------------------------------

Performing differential gene expression analysis
------------------------------------------------
