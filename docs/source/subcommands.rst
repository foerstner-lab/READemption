TRALP's subcommands
===================

align
-----

`align` includes the clipping and size filtering of the read, as well
as the actual aligning to the reference sequences.  It also generates
statistics about the steps (e.g. number of aligned reads, number of
mappings). As the result of this steps are used by the other
subcommands it has to be run before any other.

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
