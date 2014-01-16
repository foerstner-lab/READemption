TRAPL - The RNA-Seq Analysis PipeLine 
*************************************
Table of content
================

.. toctree::
   :maxdepth: 1

   index
   installation
   subcommands      
   example_analysis
   license
   versions	      

TRAPL in a nutshell
===================

*The RNA-Seq Analysis PipeLine* (*TRAPL*) is - as the name implies - a
pipeline for the computational evaluation of RNA-Seq data. It was
originally developed at the `IMIB/ZINF
<http://www.imib-wuerzburg.de/>`_ to process dRNA-Seq reads (as
introduced by Sharma *et al.*, Nature, 2010 (`Pubmed
<http://www.ncbi.nlm.nih.gov/pubmed/20164839>`_)) originating from
bacterial samples. Meanwhile is has been extended to process data
generated in different experimental setups and originating from all
domains of life and is under `active development
<https://github.com/konrad/trapl>`_. The `subcommands
<subcommands.html>`_ which are provided by command-line interface
cover read processing and aligning, coverage plot generation, gene
expression quantification as well as differential gene expression
analysis. TRAPL was applied to analyze numerous data sets. In order to
set up analyses quickly TRAPL follows the principal of *convention
over configuration*: Once the input files are copied into defined
folders no further parameters have to be given. Still, TRAPL's
behavior can be adapted to specific needs of the user. This tools is
available as open source under the `ICS <https://en.wikipedia.org/wiki/ISC_license>`_
open source license.

Download
========

TRAPL can be download from `its PyPI page
<https://pypi.python.org/pypi/trapl/>`_. Please read the
`installation instructions <installation.html>`_.

Source code
===========

The source code of TRAPL can be found at https://github.com/konrad/trapl.

Cite
====

If you apply TRALP in you data analysis please cite the following
reference: *The RNA-Seq Analysis PipeLine (TRAPL) – A tool for the
computational analysis of deep-sequencing based transcriptome data*.
Konrad U. Förstner, Jörg Vogel, Cynthia M. Sharma; (in preparation). A
`pre-preprint version <http://biorxiv.org/content/early/2014/XXX>`_ of the
manuscript is hosted at bioRxiv.

Contact
=======

For question and requests feel free to contact Konrad Förstner <konrad.foerstner@uni-wuerzburg.de>
