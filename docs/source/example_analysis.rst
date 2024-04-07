Performing example analyses
==============================

Below are two example analyses than can be performed by the user.
The first one is a single-species analysis with one species (*Salmonella* Typhimurium SL1344),
while the second one is a multi-species analysis with two species (*Staphylococcus aureus*
strain SH100 and Human Mast cells). The read files of the multi-species analysis consist of three
libraries: *Infected* has both human and staphylococcus reads, while *uninfected* has only humand and
*steady_state* only stapyhlococcus reads. Both analysis follow a similar workflow.
The main differences occur when creating a new project and running deseq,
since some information about the libraries and
the species they contain has to be provided (See details bwlow, when running the subcommands).

Single-species analysis
-----------------------

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

**Generating a project**

At first we have to create the analysis folder and its input subfolder.
All output folders, except for the *align* output folder, will be created when the corresponding subcommand is called.
For this we use the ``create`` subcommand.
Please not that all species of a project - even if it is only one species - need to be named via the *--species* flag ::

  $ reademption create --project_path READemption_analysis --species salmonella="Salmonella Typhimurium"
  Created folder "READemption_analysis" and required subfolders.
  Please copy read files into folder "READemption_analysis/input/reads" and reference sequences files into folder/s "READemption_analysis/input/salmonella_reference_sequences".

This will result in a folder structure as shown here:
::

  READemption_analysis
  ├── config.json
  ├── input
  │   ├── reads
  │   ├── salmonella_annotations
  │   └── salmonella_reference_sequences
  └── output
      └── align
          ├── alignments
          ├── index
          ├── processed_reads
          ├── reports_and_stats
          │   ├── stats_data_json
          │   └── version_log.txt
          └── unaligned_reads


**Retrieving the input data**

We have to download the reference sequence (FASTA format) as well as
the annotation file (GFF3 format) for *Salmonella* from NCBI. As we
will use the URL of *Salmonella* Typhimurium SL1344's source FTP
folder it several times we store it in an environment variable called
``FTP_SOURCE``.  

::

  $ FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/

We download the reference sequence (the chromosome and three plasmids)
in FASTA format and store them in the ``reference_sequences``
folder. The files are saved with a different suffix (``.fa`` instead
of ``.fna``) as some genome browser (e.g. IGB) will not accept them as
FASTA files otherwise.

::
   
   $ wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
   $ wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
   $ wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
   $ wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna

We have to modify the header of the FASTA files as the sequence IDs
have to be the same as the ones in the first column of the GGF3 files
(see below) to be used in the gene quantification. This will be also
necessary if both, FASTA and GFF3 files, will be loaded in the IGB.

::

   $ sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
   $ sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
   $ sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
   $ sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa

**Then we download the GFF3 files that contain the annotations**
::

   $ wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz
   $ gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz

Finally, we need the reads of the RNA-Seq libraries. To save some time
for running this examples we will work with subsampled libraries of 1M
reads each. This will the limit informative value of the results which
is acceptable as we just want to understand the workflow of the
READemption. Please be aware that READemption can perform only basic
quality trimming and adapter clipping. If this is not sufficient you
can use the `FASTX toolkit <http://hannonlab.cshl.edu/fastx_toolkit/>`_,
`cutadapt <https://code.google.com/p/cutadapt/>`_ or other tools for
the preprocessing.

::

   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
   $ wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2

We have now all the necessary data available. The input folder should
look like this now:

::

   $ ls READemption_analysis/input/* 
   READemption_analysis/input/salmonella_annotations:
   NC_016810.gff  NC_017718.gff  NC_017719.gff  NC_017720.gff
   
   READemption_analysis/input/reads:
   InSPI2_R1.fa.bz2  InSPI2_R2.fa.bz2  LSP_R1.fa.bz2  LSP_R2.fa.bz2
 
   READemption_analysis/input/salmonella_reference_sequences:
   NC_016810.fa  NC_017718.fa  NC_017719.fa  NC_017720.fa

**Processing and aligning the reads**

The first step it the read processing and mapping. Via parameters we
tell READemption to use 4 CPU (``-p 4``) and perform a poly-A-clipping
(``--poly_a_clipping``) before the mapping.

::

   $ reademption align -p 4 --poly_a_clipping --fastq --progress --project_path READemption_analysis

Once this the mapping is done the file ``read_alignment_stats.csv`` is
created which can be found in
``READemption_analysis/output/align/reports_and_stats/``. It contains
several mapping statistics for example how many reads are successfully
aligned in total and how many were aligned to each replicon. We see
that more than 98 % of the reads are mapped for each library. Sorted
and indexed alignements in BAM format are stored in
``READemption_analysis/output/align/alignments``. We could load them
into a genome browser but instead we continue with the next step.


**Generating coverage files**

In order to generate strand specific coverage files with different
normalizations we use the subcommand ``coverage``.

::

   $ reademption coverage -p 4 --project_path READemption_analysis

The sets are stored in subfolder of
``READemption_analysis/output/salmonella_coverage-raw/``, ``READemption_analysis/output/salmonella_coverage-tnoar_mil_normalized/`` and ``READemption_analysis/output/salmonella_coverage-tnoar_min_normalized/``. The most oftenly used set
is stored in ``coverage-tnoar_min_normalized``. Here the coverage
values are normalized by the total number of aligned reads (TNOAR) of
the individual library and then multiplied by the lowest TNOAR value
of all libraries.
These files could be inspected for differential
RNA-Seq (dRNA-Seq - comparing libraries with and without Terminator
Exonuclease treatment) data in order to determine transcriptional
start sites. They can be loaded in a common genome browsers like `IGB
<http://bioviz.org/igb/>`_ or `IGV
<https://www.broadinstitute.org/igv/home>`_. Keep in mind that the
coverages of the reverse strand have negative values so you have to
adapt the scaling in some genome browsers.

**Performing gene wise quantification**

In this step we want to quantify the number of reads overlapping with
the locations of the annotation entries. With the ``--features``
parameter we configure ``reademption`` to just quantify CDS, tRNA and
rRNA entries.

::

   $ reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis

After the quantification we find tables that contain the combined
counting for all entries in
``READemption_analysis/output/salmonella_gene_quanti_combined``. The
countings for mappings in sense and anti-sense are separately
listed. Besides the raw countings there are also tables for
countings normalized by the total number of reads, RPKM values and TPM (transcripts per million).


**Performing differential gene expression analysis**

To compare the gene expression of different conditions we apply the
subcommand ``deseq`` which makes use of the R library `DESeq2
<http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_. 

::

   $ reademption deseq -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 -c InSPI2,InSPI2,LSP,LSP -r 1,2,1,2 --libs_by_species salmonella=InSPI2_R1,InSPI2_R2,LSP_R1,LSP_R2 --project_path READemption_analysis

We have to tell READemption which libraries are replicates of which
condition. This is done by the parameter ``-l``, ``-c`` and ``-r`` . ``-l``
should hold a comma separated list of the libraries, ``-c`` the
corresponding conditions and ``-r`` the corresponding replicate number. In our case we have 4 libraries
(``InSPI2_R1.fa.bz2``, ``InSPI2_R2.fa.bz2``, ``LSP_R1.fa.bz2``,
``LSP_R2.fa.bz2``) and two conditions (which we call ``InSPI2`` and
``LSP``) and two times two replicates (R1 and R2 for each condition). Just to make this association easier to understand:

::
   
    libs      InSPI2_R1.fa.bz2  InSPI2_R2.fa.bz2  LSP_R1.fa.bz2  LSP_R2.fa.bz2
                 |                 |               |              |
    conds      InSPI2            InSPI2            LSP            LSP
                 |                 |               |              |
    reps         1                 2               1              2
When you call ``deseq`` it will compare all conditions with each other
and you can pick the comparison that you need. The raw ``DESeq2``
results are enriched with the original annotation information and are
stored in
``READemption_analysis/output/salmonella_deseq/deseq_with_annotations``

**Create plots**

Finally we generate plots that visualize the results of the different
steps. ``viz_align`` creates histograms of the read length
distribution for the untreated and treated reads (saved in
``READemption_analysis/output/read_lengths_viz_align/``).
It also creates an overview of how many reads map to each species
and how many reads are species cross-mapped per library (saved in
``READemption_analysis/output/all_species_viz_align/``. However, this folder can be neglected in a single species analysis).


::
   
   $ reademption viz_align --project_path READemption_analysis

``viz_gene_quanti`` visualizes the gene wise countings. In our example
you will see that - as expected - the replicates are more similar to
each other than to the libs of the other condition. It also generates
bar plots that show the distribution of reads inside the different RNA
classes.

::

    $ reademption viz_gene_quanti --project_path READemption_analysis

``viz_deseq`` generates MA-plots as well as volcano plots.


::

   $ reademption viz_deseq --project_path READemption_analysis


Multi-species analysis
----------------------


Here you will be guided trough a small example Dual RNA-seq analysis using a
publicly available RNA-Seq from the European Nucleotide Archive (ENA) that was part of a
publication by `Goldmann et
al. <https://pubmed.ncbi.nlm.nih.gov/35321877/>`_. This is a
transcriptome analysis of *Staphylococcus aureus*
strain SH100 and Human Mast cells in different
conditions. The complete analysis is publicly available at `Publisso <https://repository.publisso.de/resource/frl:6427216>`_.
Note that we use only three of the five conditions (9 instead of all 15 libraries) to make the analysis less complicated.
We will generate several output files in different
formats. The CSV (tabular separated plain text files) files can be
opened with any spreadsheet program like `LibreOffice
<https://www.libreoffice.org/>`_ or Excel. For inspecting the mappings
(in BAM format) and coverage files (wiggle format) you can use a
genome browser for example `IGB <http://bioviz.org/igb/>`_ or `IGV
<https://www.broadinstitute.org/igv/home>`_.

**Generating a project**

At first we have to create the analysis folder and its input subfolder.
All output folders, except for the *align* output folder, will be created when the corresponding subcommand is called.
For this we use the ``create`` subcommand.
Please not that all species of a project - in this case two species - need to be named via the *--species* flag ::

  $ reademption create --project_path READemption_analysis --species human="Homo sapiens" staphylococcus="Staphylococcus aureus"
  Created folder "READemption_analysis" and required subfolders.
  Please copy read files into folder "READemption_analysis/input/reads" and reference sequences files into folder/s "READemption_analysis/input/human_reference_sequences", "READemption_analysis/input/staphylococcus_reference_sequences".
This will result in a folder structure as shown here:
::

  READemption_analysis
  ├── config.json
  ├── input
  │   ├── human_annotations
  │   ├── human_reference_sequences
  │   ├── reads
  │   ├── staphylococcus_annotations
  │   └── staphylococcus_reference_sequences
  └── output
      └── align
          ├── alignments
          ├── index
          ├── processed_reads
          ├── reports_and_stats
          │   ├── stats_data_json
          │   └── version_log.txt
          └── unaligned_reads


**Retrieving the input data**

We have to download the reference sequences (FASTA format) as well as
the annotation files (GFF3 format) for both species.

|

We download the *Staphylococcus* genome to the corresponding folder and unpack it.

::

  $ wget -O READemption_analysis/input/staphylococcus_reference_sequences/staphylococcus_genome.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
  $ gunzip READemption_analysis/input/staphylococcus_reference_sequences/staphylococcus_genome.fa.gz



We download the *Staphylococcus* annotation to the corresponding folder and unpack it.

::

  $ wget -O READemption_analysis/input/staphylococcus_annotations/staphylococcus_annotation.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz
  $ gunzip READemption_analysis/input/staphylococcus_annotations/staphylococcus_annotation.gff.gz

We download the Human genome to the corresponding folder and unpack it.

::

  $ wget -O READemption_analysis/input/human_reference_sequences/human_genome.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz
  $ gunzip READemption_analysis/input/human_reference_sequences/human_genome.fa.gz

We download the Human annotation to the corresponding folder and unpack it.

::

  $ wget -O READemption_analysis/input/human_annotations/human_annotation.gff.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gff3.gz
  $ gunzip READemption_analysis/input/human_annotations/human_annotation.gff.gz

The reference *Staphylococcus sequence was saved with a different suffix (``.fa`` instead
of ``.fna``) as some genome browser (e.g. IGB) will not accept them as
FASTA files otherwise.


Finally, we need the reads of the RNA-Seq libraries. To save some time
for running this examples we will work with subsampled libraries of 10000
reads each. This will the limit informative value of the results which
is acceptable as we just want to understand the workflow of the
READemption. Please be aware that READemption can perform only basic
quality trimming and adapter clipping. If this is not sufficient you
can use the `FASTX toolkit <http://hannonlab.cshl.edu/fastx_toolkit/>`_,
`cutadapt <https://code.google.com/p/cutadapt/>`_ or other tools for
the preprocessing.

::

  $ wget https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Infected_replicate_1.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Infected_replicate_2.fq \         https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Infected_replicate_3.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Steady_state_replicate_1.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Steady_state_replicate_2.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Steady_state_replicate_3.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Uninfected_replicate_1.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Uninfected_replicate_2.fq https://raw.githubusercontent.com/Tillsa/Tillsa-2022-06-15-READemption_tutorial_data/main/Uninfected_replicate_3.fq -P READemption_analysis/input/reads

We have now all the necessary data available. The input folder should
look like this now:

::

    $ ls READemption_analysis/input/*
    input/human_annotations:
    human_annotation.gff

    input/human_reference_sequences:
    human_genome.fa

    input/reads:
    Infected_replicate_1.fq  Infected_replicate_3.fq      Steady_state_replicate_2.fq  Uninfected_replicate_1.fq  Uninfected_replicate_3.fq
    Infected_replicate_2.fq  Steady_state_replicate_1.fq  Steady_state_replicate_3.fq  Uninfected_replicate_2.fq

    input/staphylococcus_annotations:
    staphylococcus_annotation.gff

    input/staphylococcus_reference_sequences:
    staphylococcus_genome.fa


**Processing and aligning the reads**

The first step it the read processing and mapping. Via parameters we
tell READemption to use 4 CPU (``-p 4``) and perform a poly-A-clipping
(``--poly_a_clipping``) before the mapping.

::

   $ reademption align -p 4 --poly_a_clipping --project_path READemption_analysis

Once this the mapping is done the file ``read_alignment_stats.csv`` is
created which can be found in
``READemption_analysis/output/align/reports_and_stats/``. It contains
several mapping statistics for example how many reads are successfully
aligned in total and how many were aligned to each species as well as the species cross aligned reads. Sorted
and indexed alignements in BAM format are stored in
``READemption_analysis/output/align/alignments``. We could load them
into a genome browser but instead we continue with the next step.


**Generating coverage files**

In order to generate strand specific coverage files with different
normalizations we use the subcommand ``coverage``.

::

   $ reademption coverage -p 4 --project_path READemption_analysis

The sets are stored in subfolder of
``READemption_analysis/output/staphylococcus_coverage-raw/``, ``READemption_analysis/output/staphylococcus_coverage-tnoar_mil_normalized/``, ``READemption_analysis/output/staphylococcus_coverage-tnoar_min_normalized/``,
``READemption_analysis/output/human_coverage-raw/``, ``READemption_analysis/output/human_coverage-tnoar_mil_normalized/`` and ``READemption_analysis/output/human_coverage-tnoar_min_normalized/``.
The most oftenly used set is stored in ``coverage-tnoar_min_normalized``.
Here the coverage values are normalized by the total number of aligned reads (TNOAR) of
the individual library and then multiplied by the lowest TNOAR value
of all libraries.
These files could be inspected for differential
RNA-Seq (dRNA-Seq - comparing libraries with and without Terminator
Exonuclease treatment) data in order to determine transcriptional
start sites. They can be loaded in a common genome browsers like `IGB
<http://bioviz.org/igb/>`_ or `IGV
<https://www.broadinstitute.org/igv/home>`_. Keep in mind that the
coverages of the reverse strand have negative values so you have to
adapt the scaling in some genome browsers.

**Performing gene wise quantification**

In this step we want to quantify the number of reads overlapping with
the locations of the annotation entries. With the ``--features``
parameter we configure ``reademption`` to just quantify *gene* entries to save some time.

::

   $ reademption gene_quanti -p 4 --features gene --project_path READemption_analysis

After the quantification we find tables that contain the combined
counting for all entries in
``READemption_analysis/output/staphylococcus_gene_quanti_combined`` and ``READemption_analysis/output/human_gene_quanti_combined``. The
countings for mappings in sense and anti-sense are separately
listed. Besides the raw countings there are also tables for
countings normalized by the total number of reads, RPKM values and TPM (transcripts per million).


**Performing differential gene expression analysis**

To compare the gene expression of different conditions we apply the
subcommand ``deseq`` which makes use of the R library `DESeq2
<http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_.

::

   $ reademption deseq -l Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Steady_state_replicate_1,Steady_state_replicate_2,Steady_state_replicate_3,Uninfected_replicate_1,Uninfected_replicate_2,Uninfected_replicate_3 -c infected,infected,infected,steady_state,steady_state,steady_state,uninfected,uninfected,uninfected -r 1,2,3,1,2,3,1,2,3 --libs_by_species human="Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Uninfected_replicate_1,Uninfected_replicate_2,Uninfected_replicate_3" staphylococcus="Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Steady_state_replicate_1,Steady_state_replicate_2,Steady_state_replicate_3" --size_factor=species --project_path READemption_analysis

We have to tell READemption which libraries are replicates of which
condition. This is done by the parameter ``-l``, ``-c`` and ``-r`` . ``-l``
should hold a comma separated list of the libraries, ``-c`` the
corresponding conditions and ``-r`` the corresponding replicate number.
In our case we have 9 libraries (``Infected_replicate_1``, ``Infected_replicate_2``, ``Infected_replicate_3``, ``Steady_state_replicate_1``, ``Steady_state_replicate_2``, ``Steady_state_replicate_3``, ``Uninfected_replicate_1``, ``Uninfected_replicate_2``, ``Uninfected_replicate_3``)
and three conditions (which we call ``infected``, ``steady_state`` and ``uninfected``) and three times three replicates (R1, R2 and R3 for each condition). Just to make this association easier to understand:

::

    libs      Infected_replicate_1    Infected_replicate_2    Infected_replicate_3    Steady_state_replicate_1    Steady_state_replicate_2    Steady_state_replicate_3    Uninfected_replicate_1    Uninfected_replicate_2    Uninfected_replicate_3
                 |                              |                       |                        |                        |                               |                         |                         |                         |
    conds     infected                      infected                infected               steady_state             steady_state                    steady_state                uninfected                uninfected                uninfected
                 |                              |                       |                        |                        |                               |                         |                         |                         |
    reps         1                              2                       3                        1                        2                               3                         1                         2                         3
Because we set the --size factor to species and set the species for each lib via --libs_by_species,
when you call ``deseq`` it will calculate the size factors for normalization based on the reads of the current species.
Deseq will compare all the possible combinations of the libraries of a species
and you can pick the comparison that you need. The raw ``DESeq2``
results are enriched with the original annotation information and are
stored in
``READemption_analysis/output/staphyloccus_deseq/deseq_with_annotations`` and ``READemption_analysis/output/human_deseq/deseq_with_annotations``

**Create plots**

Finally we generate plots that visualize the results of the different
steps. ``viz_align`` creates histograms of the read length
distribution for the untreated and treated reads (saved in
``READemption_analysis/output/read_lengths_viz_align/``).
It also creates an overview of how many reads map to each species and how many reads are species cross-mapped per library.

::

   $ reademption viz_align --project_path READemption_analysis

``viz_gene_quanti`` visualizes the gene wise countings. In our example
you will see that - as expected - the replicates are more similar to
each other than to the libs of the other condition. It also generates
bar plots that show the distribution of reads inside the different RNA
classes.

::

   $ reademption viz_gene_quanti --project_path READemption_analysis

``viz_deseq`` generates MA-plots as well as volcano plots.

::

   $ reademption viz_deseq --project_path READemption_analysis
