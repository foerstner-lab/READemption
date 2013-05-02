RAPL - The RNA-Seq Analysis PipeLine 
====================================

.. toctree::
   :maxdepth: 2

What is RAPL?
=============

*RNA-Seq Analysis PipeLine* (*RAPL*) is - as the name implies - a pipeline
for the computational analysis of RNA-Seq data. It was originally
developed to process dRNA-Seq reads (see Sharma *et al.*, Nature,
2010) originating from bacterial samples. Meanwhile is was extended to
process data generated in different experimental setups and
originating from all domain of life. It was applied to analyze
numerous data sets.

RAPL follows the principal of *convention over configuration*

RALP's subcommands
==================

align
-----

`align` includes the clipping and size filtering of the read, as well
as the actual aligning to the reference sequences.  It also generates
statistic about the steps (e.g. number of aligned reads, number of
mappings). The result of this steps are used by the other subcommands
so it has to be run before any other.

coverage
--------

`coverage` generates strand specific coverage files in wiggle format
based on the read alignments. The wiggle files can be viewed in common
genome browser like the `Integrated genome browser (IGB)
<http://bioviz.org/>`_ or the `Integrative genome viewer (IGV)
<https://www.broadinstitute.org/software/igv/>`_.

gene_quanti
-----------

With `gene_quanti` the number of reads to each gene is counted and the
results are combined in tables.

deseq
-----

Differential gene expression can be performed using `deseq` which will
run a `DESeq <http://www-huber.embl.de/users/anders/DESeq/>`_ analyses for all possible combinations.

viz_align
---------

`viz_align` generates plots that visualize selected features from the
alignment results.

viz_gene_quanti
---------------

viz_deseq
---------

Installation
============

Requirements
------------

RAPL was developed using Python 3.3 and for best performance the user
is advised to use this version, too. Any other Python 3 version should
work as well. Also Python 2.7 can be used if the library XXX is
installed. In any case, the third party module `pysam
<https://code.google.com/p/pysam>`_ is required. The short read mapper


`segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_

Matplotlib for `viz_align`, `viz_gene_quanti`, `viz_deseq`

R and DESeq for differential gene expression analysis

Installing in a fresh Ubuntu environment
----------------------------------------

Amazon AWS, Ubuntu live system

Global installation
-------------------

Installation in the home directory of the user
----------------------------------------------

Installation in a venv
----------------------

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

Cite
====

If you use RALP for your research please cite the following reference:
*RNA-Seq Analysis PipeLine (RAPL) – A tool for the computational
analysis of deep-sequencing based transcriptome data;
Konrad U. Förstner, Jörg Vogel, Cynthia M. Sharma; (in preparation)

Source code
===========

The source code of RAPL can be found at https://github.com/konrad/rapl.

License
=======

RAPL is open source software and available under the ISC license.

Copyright (c) 2013, Konrad Förstner <konrad.foerstner@uni-wuerzburg.de>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

Versions/Change log
===================
