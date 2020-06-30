READemption's subcommands
=========================

In general the subcommands need at least one argument - the analysis
folder.

create
------

``create`` generates the required folder structure for input and
output files. Once these folders are created the input files have to
be placed into the correct locations. As a minimal requirement,
RNA-Seqs reads in FASTA format (can be compressed with ``bzip2`` or
``gzip``) must be placed in ``input/reads`` and the reference sequence
in FASTA format must be copied or linked in
``input/reference_sequences``. For the command ``gene_quanti``
annotation files in GFF3 format have to be put into
``input/annotations``.

.. argparse::
   :filename: ../bin/reademption
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
several mapping statistics. The folder
``output/align/reports_and_stats/stats_data_json/`` contains files
with the original countings in JSON format. Please be aware that
READemption can perform only basic quality trimming and adapter
clipping. If this is not sufficient you can use the `FASTX toolkit
<http://hannonlab.cshl.edu/fastx_toolkit/>`_, `cutadapt
<https://code.google.com/p/cutadapt/>`_ or other tools for the
preprocessing.

.. argparse::
   :filename: ../bin/reademption
   :prog: reademption
   :func: create_parser
   :path: align

coverage
--------

`coverage` generates strand specific coverage files in `wiggle format
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
   :filename: ../bin/reademption
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
raw read countings, Transcripts per Million (TPM) normalized read counts
as well as Reads per Kilobase Million (RPKM) normalized read counts.
The results of the different counting methods are provided in separate CSV-tables.

.. argparse::
   :filename: ../bin/reademption
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
   :filename: ../bin/reademption
   :prog: reademption
   :func: create_parser
   :path: deseq


viz_align
---------

``viz_align`` plots histograms of the read length distributions of the
reads before and after the read clipping.

.. argparse::
   :filename: ../bin/reademption
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
   :filename: ../bin/reademption
   :prog: reademption
   :func: create_parser
   :path: viz_gene_quanti


viz_deseq
---------

``viz_deseq`` generates MA-plots of the comparison (log2 fold changes
vs. the base mean) as well as volcano plots (log2 fold changes
vs. p-values / adjusted p-values).

.. argparse::
   :filename: ../bin/reademption
   :prog: reademption
   :func: create_parser
   :path: viz_deseq
