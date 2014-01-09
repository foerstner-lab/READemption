= About =

The RNA-Seq Analysis PipeLine (TRAPL) is - as the name implies - a
pipeline for the computational evaluation of RNA-Seq data. It was
originally developed at the IMIB/ZINF to process dRNA-Seq reads (as
introduced by Sharma et al., Nature, 2010 originating from bacterial
samples. Meanwhile is has been extended to process data generated in
different experimental setups and originating from all domains of life
and is under active developtment. The subcommands which are provided
by command-line interface cover read processing and aligning, coverage
plot generation, gene expression quantification as well as
differential gene expression analysis. TRAPL was applied to analyze
numerous data sets. In order to set up analyses quickly TRAPL follows
the principal of "convention over configuration": Once the input files
are copied into defined folders no further parameters have to be
given. Still, TRAPL's behavior can be adapted to specific needs of the
user.

= License =

ICS - see LICENSE.txt 

== Development ==

- If possible follow the principal of "convention over
  configuration". This means input file are into a fixed location and
  the result file are placed in fixed location.

- The classes should be path agnostic as far a possible. The controller
  is taking care of that and call them adequately.
