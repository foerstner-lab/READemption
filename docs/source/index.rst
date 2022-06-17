READemption - A RNA-Seq Analysis Pipeline 
*****************************************
Table of content
================

.. toctree::
   :maxdepth: 1

   installation
   live_and_installation_image
   docker_image
   subcommands
   fragment_building
   example_analysis
   troubleshooting
   license
   versions	      

READemption in a nutshell
=========================

READemption is a pipeline for the computational evaluation of
RNA-Seq data. It was originally developed to process dRNA-Seq reads
(as introduced by Sharma et al., Nature, 2010 (`Pubmed
<http://www.ncbi.nlm.nih.gov/pubmed/20164839>`_)) originating from
bacterial samples. Meanwhile it has been extended to process data
generated in different experimental setups and from all domains of
life.
READemption features handling of:

- **single and multi-species projects**

- **single or paired-end reads**

The `functions <subcommands.html>`_ which are accessible via a
command-line interface cover:

- **read processing**

- **aligning**

- **coverage calculation**

- **gene expression quantification**

- **differential gene expression analysis**

- **visualization**

In order to set up and
perform analyses quickly READemption follows the principal of
*convention over configuration*: Once the input files are
copied/linked into defined folders no further parameters have to be
given. Still, READemption's behavior can be adapted to specific needs
of the user by parameters.

Download
========

READemption can be download from `its PyPI page
<https://pypi.python.org/pypi/READemption/>`_. Please read the
`installation instructions <installation.html>`_.

Source code
===========

The source code of READemption can be found at https://github.com/foerstner-lab/READemption.

Cite
====

If you apply READemption in your data analysis please cite the
following publication: *READemption – A tool for the computational
analysis of deep-sequencing-based transcriptome data*.
Konrad U. Förstner, Jörg Vogel, Cynthia M. Sharma. 2014,
Bioinformatics, Aug 13. `Fulltext
<http://bioinformatics.oxfordjournals.org/content/30/23/3421>`_,
`Pre-print at bioRxiv
<http://biorxiv.org/content/early/2014/05/19/003723>`_.

Contact
=======

For question and requests feel free to contact `Konrad Förstner
<http://konrad.foerstner.org/>`_ <foerstner@zbmed.de>
