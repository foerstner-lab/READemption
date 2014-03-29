READemption's subcommands
=========================

create
------

`create` generates the required folder structure for input and output
files which looks likes this:

::

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
    │   │   └── used_rapl_version.txt
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

align
-----

`align` performs the clipping and size filtering of the reads, as well
as the actual aligning to the reference sequences. It also generates
statistics about the steps (e.g. number of aligned reads, number of
mappings). As the result of this steps are needed by the other
subcommands it has to be run before any other. It requires reads in
FASTA format (or counterparts compressed with ``gzip`` or ``bzip2``)
and reference sequences in FASTA format (one sequence per
file). `align` generates the read alignment in BAM format (`*.bam`)
and also index files for those (`*.bai`). Is also stores unmapped

=> report
::

 positional arguments:
   project_path          Path of the project folder. If none is given the
                         current directory is used.

 optional arguments:
  -h, --help            show this help message and exit
  --min_read_length MIN_READ_LENGTH, -l MIN_READ_LENGTH
                        Minimal read length after clipping
  --processes PROCESSES, -p PROCESSES
                        Number of processes that should be used.
  --segemehl_accuracy SEGEMEHL_ACCURACY, -a SEGEMEHL_ACCURACY
                        Segemehl's minimal accuracy (in %) (default 95).
  --segemehl_evalue SEGEMEHL_EVALUE, -e SEGEMEHL_EVALUE
                        Segemehl's maximal e-value. (default 5.0)
  --segemehl_bin SEGEMEHL_BIN, -s SEGEMEHL_BIN
                        Segemehl's binary path.
  --split, -S           Run segemehl with read splitting
  --poly_a_clipping, -c
                        Perform polyA tail clipping. This option cannot be
                        used for paired-end reads.
  --force, -f           Overwrite existing files.
  --progress, -P        Show progress of the segemehl mapping.
  --paired_end, -r      Use this if reads are originating from a paired-end
                        sequencing. The members of a pair must be marked with
                        '_p1' and '_p2' in front of the file type suffixes
                        (e.g. 'my_sample_p1.fa' and 'my_sample_p2.fa' or
                        'my_sample_p1.fa.bz2' and 'my_sample_p2.fa.bz2'). This
                        option cannot be use with polyA tail clipping.

coverage
--------

`coverage` generates strand specific coverage files in `wiggle format
<http://genome.ucsc.edu/goldenPath/help/wiggle.html>`_ based on the
read alignments. These wiggle files can be viewed in common genome
browser like the `Integrated genome browser (IGB)
<http://bioviz.org/>`_ or the `Integrative genome viewer (IGV)
<https://www.broadinstitute.org/software/igv/>`_. Three sets of wiggle
files will be generated: raw couting values without normalization
(located in the folder `coverage-raw`), normalized by the total number
of aligned reads (abriviated as tnoar) and the multiplied by the
lowest number of aligned reads of all considered libraries (in folder
`coverage-tnoar_min_normalized`) as well as normalized by the total
number of aligned reads and multiplied by one million
(`coverage-tnoar_mil_normalized`). The different normalisations make a
visual semi-quantitative comparative possible and enable to perform
transcription start site analysis (e.g. using tools like `TSSAR
<http://nylon.tbi.univie.ac.at/TSSAR/>`_). For each library and set
there will be coverage files for the forward and the reverse
strand. The coverages for the forward strand have positive while the
one for the reverse stand have negative values in order to make a
visual discrimanation easy. Per default all reads and each position of
them will be considered. To calculate the coverages only based on
uniquely aligned read use the ``--unique_only`` parameter. If only the
first base should be considered add ``--first_base_only``. Reads are
aligned to multiple location will account only in fraction to the
values of the different positions. For example a read that is mapped
to three different location will contribute a value of 1/3 to each of
the nucleotiedes of these positions. To turn off this behavior use
``--skip_read_count_splitting``.

:: 

 positional arguments:
  project_path          Path of the project folder. If none is given the
                        current directory is used.

 optional arguments:
  -h, --help            show this help message and exit
  --unique_only, -u     Use uniquely aligned reads only.
  --processes PROCESSES, -p PROCESSES
                        Number of processes that should be used.
  --skip_read_count_splitting, -s
                        Do not split the read counting between different
                        alignings. Default is to do the splitting.
  --first_base_only, -b
                        Only the first bases 5' base of each read aligning is
                        taken into account.
  --force, -f           Overwrite existing files.

gene_quanti
-----------

With `gene_quanti` the number of reads to each gene is counted and the
results are combined in tables. 

- IDs must be the same
::

 positional arguments:
  project_path          Path of the project folder. If none is given the
                        current directory is used.

 optional arguments:
  -h, --help            show this help message and exit
  --min_overlap MIN_OVERLAP, -o MIN_OVERLAP
                        Minimal read-annotation-overlap (in nt) (default 1)
  --skip_norm_by_alignment_freq
  --skip_norm_by_overlap_freq
  --skip_antisense, -a
  --processes PROCESSES, -p PROCESSES
                        Number of processes that should be used.
  --features ALLOWED_FEATURES, -t ALLOWED_FEATURES
                        Comma separated list of features that should be
                        considered (e.g. gene, cds, region, exon). Other
                        feature will be skipped. If not specified all features
                        will be considered.
  --unique_only, -u     Use uniquely aligned reads only.
  --pseudocounts, -c    Add a pseudocount of 1 to each gene.
  --force, -f           Overwrite existing files.

deseq
-----

Differential gene expression can be performed using `deseq` which will
run a `DESeq <http://www-huber.embl.de/users/anders/DESeq/>`_ analyses for all possible combinations.

::

 positional arguments:
  project_path          Path of the project folder. If none is given the
                        current directory is used.

 optional arguments:
  -h, --help            show this help message and exit
  --libs LIBS, -l LIBS  Comma separated list of libraries.
  --conditions CONDITIONS, -c CONDITIONS
                        Comma separated list of condition in the same order as
                        their corresponding libraries.
  --no_replicates, -r


viz_align
---------

`viz_align` plots histograms of the read length distributions of the
reads before and after the read clipping.

viz_gene_quanti
---------------

`viz_gene_quanti` creates scatterplots in with the raw gene wise
quantification values are compared for each library pair
(all-against-all). For each comparison the `pearson correllation
<https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient>`_
(`r`) coefficiant is.

viz_deseq
---------
