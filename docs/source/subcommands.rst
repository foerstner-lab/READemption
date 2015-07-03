READemption's subcommands
=========================

In general the subcommands need at least one argument - the analysis
folder. If this is not given READemption assumes that the current
folder is the analysis folder.

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

::

   usage: reademption create [-h] project_path

   positional arguments:
     project_path  Name/path of the project.

   optional arguments:
     -h, --help    show this help message and exit

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

::

  usage: reademption align [-h] [--min_read_length MIN_READ_LENGTH]
                           [--processes PROCESSES]
                           [--segemehl_accuracy SEGEMEHL_ACCURACY]
                           [--segemehl_evalue SEGEMEHL_EVALUE]
                           [--segemehl_bin SEGEMEHL_BIN] [--paired_end]
                           [--split] [--poly_a_clipping] [--realign]
                           [--keep_original_alignments] [--lack_bin LACK_BIN]
                           [--fastq] [--check_for_existing_files] [--progress]
                           [--crossalign_cleaning CROSSALIGN_CLEANING_STRING]
                           [project_path]
  
  positional arguments:
    project_path          Path of the project folder. If none is given the
                          current directory is used.
  
  optional arguments:
    -h, --help            show this help message and exit
    --min_read_length MIN_READ_LENGTH, -l MIN_READ_LENGTH
                          Minimal read length after clipping (default 12).
                          Should be higher for eukaryotic species.
    --processes PROCESSES, -p PROCESSES
                          Number of processes that should be used (default 1).
    --segemehl_accuracy SEGEMEHL_ACCURACY, -a SEGEMEHL_ACCURACY
                          Segemehl's minimal accuracy (in %) (default 95).
    --segemehl_evalue SEGEMEHL_EVALUE, -e SEGEMEHL_EVALUE
                          Segemehl's maximal e-value (default 5.0).
    --segemehl_bin SEGEMEHL_BIN, -s SEGEMEHL_BIN
                          Segemehl's binary path (default 'segemehl.x').
    --paired_end, -P      Use this if reads are originating from a paired-end
                          sequencing. The members of a pair must be marked with
                          '_p1' and '_p2' in front of the file type suffixes
                          (e.g. 'my_sample_p1.fa' and 'my_sample_p2.fa' or
                          'my_sample_p1.fa.bz2' and 'my_sample_p2.fa.bz2'). This
                          option cannot be use with polyA tail clipping.
    --split, -S           Run segemehl with read splitting.
    --poly_a_clipping, -c
                          Perform polyA tail clipping. This option cannot be
                          used for paired-end reads.
    --realign, -r         Perform realignment of unmapped reads using 'lack'.
    --keep_original_alignments, -k
                          Only used with --realign/-r. Keep the alignment file
                          of the primary mapper (segemehl) and the realigner
                          (lack) after merging.
    --lack_bin LACK_BIN, -L LACK_BIN
                          Lack's binary path (default 'lack.x').
    --fastq, -q           Input reads are in FASTQ not FASTA format.
    --min_phred_score MIN_PHRED_SCORE, -Q MIN_PHRED_SCORE
                          Minimal Phred score. Works only if read are given in
                          FASTQ format. As soon as a based drop below this value
                          it and all the nucleotides downstream of it will be
                          trimmed off.
    --adapter ADAPTER, -A ADAPTER
                          Adapter sequence. If it is found in a read it and all
                          the nucleotides downstream will be trimmed off.
    --check_for_existing_files, -f
                          Check for existing files (e.g. from a interrupted
                          previous run) and do not overwrite them if they exits.
                          Attention! You have to take care that there are no
                          partially generated files left!
    --reverse_complement, -R
                          Map reverse complement of the input reads.
    --progress, -g          Show progress of the segemehl mapping.
    --crossalign_cleaning CROSSALIGN_CLEANING_STRING, -x CROSSALIGN_CLEANING_STRING
                          Remove reads that are cross-mapped to replicons of
                          different species. To associated species and replicons
                          give a string in the following format: '<ORG_NAME_1>:<
                          org_1_repl1>,<org_1_repl2>,..,<org_1_repl_n>;<ORG_NAME
                          _2>:<org_2_repl1>,<org_2_repl2>,..,<org_2_repl_n>'

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
contribute a value of 1/3 to each of the nucleotiedes of these
positions. To turn off this behavior use
``--skip_read_count_splitting``.

:: 

  usage: reademption coverage [-h] [--unique_only] [--normalize_by_uniquely]
                              [--processes PROCESSES]
                              [--skip_read_count_splitting] [--first_base_only]
                              [--check_for_existing_files]
                              [project_path]

  positional arguments:
    project_path          Path of the project folder. If none is given the
                          current directory is used.
  
  optional arguments:
    -h, --help            show this help message and exit
    --unique_only, -u     Use uniquely aligned reads only.
    --normalize_by_uniquely, -U
                          Normalize by the number of uniquely aligned reads. By
                          default the normalization is done based on the total
                          number of aligned reads even if only uniquely aligned
                          reads are used for the coverage calculation.
    --processes PROCESSES, -p PROCESSES
                          Number of processes that should be used (default 1).
    --skip_read_count_splitting, -s
                          Do not split the read counting between different
                          alignings. Default is to do the splitting.
    --non_strand_specific, -d
                          Do not distict between the coverage of the forward and
                          reverse strand but sum them to a single value for each
                          base.
    --first_base_only, -b
                          Only the first bases 5' base of each read aligning is
                          taken into account.
    --check_for_existing_files, -f
                          Check for existing files (e.g. from a interrupted
                          previous run) and do not overwrite them if they exits.
                          Attention! You have to take care that there are no
                          partially generated files left!

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
overlaps are counted and separately listed.

::

  usage: reademption gene_quanti [-h] [--min_overlap MIN_OVERLAP]
                                 [--no_count_split_by_alignment_no]
                                 [--no_count_splitting_by_gene_no]
                                 [--skip_antisense] [--processes PROCESSES]
                                 [--features ALLOWED_FEATURES] [--unique_only]
                                 [--pseudocounts] [--check_for_existing_files]
                                 [project_path]

  positional arguments:
    project_path          Path of the project folder. If none is given the
                          current directory is used.
  
  optional arguments:
    -h, --help            show this help message and exit
    --min_overlap MIN_OVERLAP, -o MIN_OVERLAP
                          Minimal read-annotation-overlap (in nt) (default 1).
    --no_count_split_by_alignment_no, -n
                          Do not split read countings by the number of
                          alignments a read has. By default this count splitting
                          is performed.
    --no_count_splitting_by_gene_no, -l
                          Do not split read countings by the number of genes it
                          overlaps with. By default this count splitting is
                          performed.
    --skip_antisense, -a  Do not count anti-sense read-gene-overlaps. By default
                          sense and anti-sense overlaps are counted and
                          separately reported.
    --non_strand_specific
                          Use countings of reads overlapping with a gene on both
                          strands and sum them up.
    --processes PROCESSES, -p PROCESSES
                          Number of processes that should be used (default 1).
    --features ALLOWED_FEATURES, -t ALLOWED_FEATURES
                          Comma separated list of features that should be
                          considered (e.g. gene, cds, region, exon). Other
                          feature will be skipped. If not specified all features
                          will be considered.
    --unique_only, -u     Use uniquely aligned reads only.
    --pseudocounts, -c    Add a pseudocount of 1 to each gene.
    --check_for_existing_files, -f
                          Check for existing files (e.g. from a interrupted
                          previous run) and do not overwrite them if they exits.
                          Attention! You have to take care that there are no
                          partially generated files left!

deseq
-----

Differential gene expression can be performed using ``deseq`` which
will run a `DESeq2 <http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html>`_
analyses for all possible combinations of conditions. To allocated the
conditions to the libraries use the ``--libs`` and ``--conditions``
parameters (e.g. ``--libs
SamA_R1.fa,SamA_R2.fa,SamB_R1.fa,SamB_R2.fa --conditions
SamA,SamA,SamB,SamB``).

::

  usage: reademption deseq [-h] --libs LIBS --conditions CONDITIONS
                           [--cooks_cutoff_off]
                           [project_path]
  
  positional arguments:
    project_path          Path of the project folder. If none is given the
                          current directory is used.
  
  optional arguments:
    -h, --help            show this help message and exit
    --libs LIBS, -l LIBS  Comma separated list of libraries.
    --conditions CONDITIONS, -c CONDITIONS
                          Comma separated list of condition in the same order as
                          their corresponding libraries.
    --cooks_cutoff_off, -k


viz_align
---------

``viz_align`` plots histograms of the read length distributions of the
reads before and after the read clipping.

::

  usage: reademption viz_align [-h] [project_path]

  positional arguments:
    project_path  Path of the project folder. If none is given the current
                  directory is used.

  optional arguments:
    -h, --help    show this help message and exit

viz_gene_quanti
---------------

``viz_gene_quanti`` creates scatterplots in which the raw gene wise
quantification values are compared for each library pair
(all-against-all). For each comparison the `pearson correllation
<https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient>`_
(`r`) coefficiant is. Additionally, bar charts that visualize the
distribution of the read counting of the different annotation classes
are plotted.

::

  usage: reademption viz_gene_quanti [-h] [project_path]

  positional arguments:
    project_path  Path of the project folder. If none is given the current
                  directory is used.

  optional arguments:
    -h, --help    show this help message and exit

viz_deseq
---------

``viz_deseq`` generates MA-plots of the comparison (log2 fold changes
vs. the base mean) as well as volcano plots (log2 fold changes
vs. p-values / adjusted p-values).

::

  usage: reademption viz_deseq [-h] [project_path]

  positional arguments:
    project_path  Path of the project folder. If none is given the current
                  directory is used.

  optional arguments:
    -h, --help    show this help message and exit
  --max_pvalue MAX_PVALUE
                          Maximum adjusted p-value for genes considered to be
                          regulated. Genes with adjusted p-values below will be
                          marked red. (default 0.05)
