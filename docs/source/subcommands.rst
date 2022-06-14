READemption's subcommands
=========================

**Results separately for each species**

In general the subcommands need at least one argument - the analysis
folder. The only exception is the ``create`` subcommand, which needs the ``species`` of a project.
After the initial project creation and the subsequent alignment,
which is (except for the statistics) species agnostic all other subcommands
(``coverage``, ``gene_quanti``, ``deseq``, ``viz_align``, ``viz_gene_quanti``, ``viz_deseq``) present
their results by species.


**Exclude or include species cross-mapped reads**

Usually a multi-species project will contain cross-mapped reads that map to different species.
These reads are excluded from countings and normalization factors of the subcommands ``coverage`` and ``gene_quanti`` per default.
Since the other subcommands (``deseq``, ``viz_gene_quanti``, ``viz_deseq``) are based on ``gene_quanti``,
species cross-mapped reads are also excluded from these subcommands.
This default behaviour can be turned off separately for countings [--count_cross_aligned_reads] and normalization [--normalize_cross_aligned_reads_included],
when executing ``coverage`` or ``gene_quanti``.

**Fragment building for paired-end reads**

When using paired-end reads, READemption automatically builds fragments from aligned reads
as part of the ``align`` subcommand.
This behaviour, if not needed, can be turned off [--no_fragment_building] to save time.
To continue skipping fragment building, the flag [--no_fragment_building] also needs to be set during
``coverage`` or ``gene_quanti``. Otherwise the fragments will be build and used by ``coverage`` and ``gene_quanti``,
even if fragments were not built during ``align``.

create
------

``create`` generates the required folder structure for input and
output files. A ``folder prefix`` and a ``display name`` must be assigned to each ``species``,
even if the project contains only a single species.
Once these folders are created the input files have to
be placed into the correct locations. A ``reference_sequences`` folder and an
``annotation`` folder will be created for each species, while the ``reads`` folder
is common for all species.

The following example creates the input folders for a project with three different species::

  $ reademption create --project_path reademption_analysis_triple \
        --species \
        human="Homo sapiens" \
        staphylococcus="Staphylococcus aureus" \
        influenza="Influenza A"

Created input folder structure:

| reademption_analysis_triple
| ├── config.json
| ├── input
| │      ├── human_annotations
| │      ├── human_reference_sequences
| │      ├── influenza_annotations
| │      ├── influenza_reference_sequences
| │      ├── reads
| │      ├── staphylococcus_annotations
| │      └── staphylococcus_reference_sequences
| └── output
|          └── align
|                   ├── alignments
|                   ├── index
|                   ├── processed_reads
|                   ├── reports_and_stats
|                   │      ├── stats_data_json
|                   │      └── version_log.txt
|                   └── unaligned_reads



As a minimal requirement,
RNA-Seqs reads in FASTA format (can be compressed with ``bzip2`` or
``gzip``) must be placed in ``input/reads`` and the reference sequences of the different species
in FASTA format must be copied or linked to their corresponding species input folder
``input/[species]_reference_sequences``. For the command ``gene_quanti``
annotation files in GFF3 format have to be put into
``input/[species]_annotations``.
Please note that a single species project will have only one annotations and one reference_sequences folder.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: create

align
-----

``align`` performs the clipping and size filtering of the reads, as
well as the actual aligning to the reference sequences. It also
generates statistics about the steps (e.g. number of aligned reads,
number of mappings). As the result of this steps are needed by the
other subcommands it has to be run before the others. It requires
reads in FASTA or FASTQ format (or counterparts compressed with
``gzip`` or ``bzip2``) and reference sequences in FASTA
format. ``align`` generates the read alignments in BAM format
(``*.bam``) and also index files for those (``*.bam.bai``). Is also
stores unmapped reads so that they can be inspected e.g. to search for
contaminations. The file
``output/align/reports_and_stats/read_alignment_stats.csv`` lists
several mapping statistics for the whole project, separated by species and per chromosome.
The folder ``output/align/reports_and_stats/stats_data_json/`` contains files
with the original countings in JSON format. Please be aware that
READemption can perform only basic quality trimming and adapter
clipping. If this is not sufficient you can use the `FASTX toolkit
<http://hannonlab.cshl.edu/fastx_toolkit/>`_, `cutadapt
<https://code.google.com/p/cutadapt/>`_ or other tools for the
preprocessing.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: align

coverage
--------

``coverage`` generates strand specific coverage files in `wiggle format
<http://genome.ucsc.edu/goldenPath/help/wiggle.html>`_ based on the
read alignments. These wiggle files can be viewed in common genome
browser like the `Integrated genome browser (IGB)
<http://bioviz.org/>`_ or the `Integrative genome viewer (IGV)
<https://www.broadinstitute.org/software/igv/>`_. Three sets of wiggle
files will be generated: raw counting values without normalization
(located in the folder `coverage-raw`), normalized by the total number
of aligned reads (abbreviated as tnoar) and the multiplied by the
lowest number of aligned reads of all considered libraries (in folder
`coverage-tnoar_min_normalized`) as well as normalized by the total
number of aligned reads and multiplied by one million
(`coverage-tnoar_mil_normalized`). The different normalizations make a
visual semi-quantitative comparative possible and enable to perform
transcription start site analysis (e.g. using tools like `TSSPredator
<http://www-ps.informatik.uni-tuebingen.de/itNew/?page_id=1860>`_). For
each library and set there will be coverage files for the forward and
the reverse strand. The coverages for the forward strand have positive
values while the one for the reverse stand have negative values in
order to make a visual discrimination easy. Per default all reads and
each position of them will be considered. To calculate the coverages
only based on uniquely aligned read use the ``--unique_only``
parameter. If only the first base should be considered add
``--first_base_only``. Reads are aligned to multiple location will
account only in fraction to the values of the different positions. For
example a read that is mapped to three different location will
contribute a value of 1/3 to each of the nucleotides of these
positions. To turn off this behavior use
``--skip_read_count_splitting``.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: coverage

gene_quanti
-----------

With ``gene_quanti`` the number of reads overlapping with each of the
annotation entries is counted and the results are combined in
tables. At least one GGF3 file with annotations has to be placed in
``input/annotations``. The sequence ID of the sequenced must be
precisely the same as the IDs used in the reference sequence FASTA
files. To specify the feature classes (the third column in the GFF3
file e.g. CDS, gene, rRNA, tRNA) that should be quantified the
parameter ``--features`` can be used. Otherwise countings for all
annotation entries are generated. Per default sense and anti-sense
overlaps are counted and separately listed. ``gene_quanti`` provides, besides the
raw read countings, Transcripts per Million (TPM) normalized read counts,
normalized by the total number of aligned reads of the given library (TNOAR),
as well as Reads per Kilobase Million (RPKM) normalized read counts
The results of the different counting methods are provided in separate CSV-tables.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: gene_quanti

deseq
-----

Differential gene expression can be performed using ``deseq`` which
will run a `DESeq2 <http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_
analyses for all possible combinations of conditions. To allocated the
conditions to the libraries use the ``--libs`` and ``--conditions``
parameters (e.g. ``--libs
SamA_R1.fa,SamA_R2.fa,SamB_R1.fa,SamB_R2.fa --conditions
SamA,SamA,SamB,SamB``).

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: deseq


viz_align
---------

``viz_align`` plots histograms of the read length distributions of the
reads before and after the read clipping.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: viz_align

viz_gene_quanti
---------------

``viz_gene_quanti`` creates scatterplots in which the raw gene wise
quantification values are compared for each library pair
(all-against-all). For each comparison the `pearson correllation
<https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient>`_
(`r`) coefficiant is. Additionally, bar charts that visualize the
distribution of the read counting of the different annotation classes
are plotted.

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: viz_gene_quanti


viz_deseq
---------

``viz_deseq`` generates MA-plots of the comparison (log2 fold changes
vs. the base mean) as well as volcano plots (log2 fold changes
vs. p-values / adjusted p-values).

.. argparse::
   :filename: reademption
   :prog: reademption
   :func: create_parser
   :path: viz_deseq
