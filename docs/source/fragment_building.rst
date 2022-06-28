Fragment building for paired-end reads
======================================

READemption builds fragments for paired-end reads.
Some sequencing protocols generate paired-end reads, where a template fragment is sequenced from
both sides, resulting in a *read one* in 5'-3' orientation and a *read two* in 3'-5' orientation.
After both reads have been aligned, a start and an end position can be calculated for the template.
READemption derives these positions from the alignments of a pair and writes them to a BAM file.
The resulting BAM file will consist of single-end reads that present the original template fragments.
If a paired-end read could not be properly mapped in pair, i.e. only one read of a pair could be mapped,
the mapped read will also be written to the BAM file, but its initial start and stop positions will be preserved.
The figures below show different possible alignment layouts for the two reads of a pair and their calculated fragment.

Fragment building for different layouts of pairs that map to the forward strand
-------------------------------------------------------------------------------

**A\) Pair in order, no overlap**


       start 1                         end 1

       -----------------read 1-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                                                      start 2                        end 2

                                                                      <---------------read 2--------------

       -----------------------------calculated fragment--------------------------------------->


|

|

|

**B\) Pair in order, overlap**


       start 1                         end 1

       -----------------read 1-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                            start 2                        end 2

                                            <---------------read 2--------------

       ------------------calculated fragment------------------------->

|

|

|

**C\) Pair in order, identical overlap**


                                            start 1                         end 1

                                            -----------------read 1-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                            start 2                        end 2

                                            <---------------read 2--------------

                                            ----calculated fragment--->

|

|

|

**D\) Pair in reverse order, both reads exceed each other (template fragment length is smaller than read length)**


                                                   start 1                         end 1

                                                   -----------------read 1-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                     start 2                        end 2

                                     <---------------read 2--------------

    calculated fragment:              ------------------------>

|

|

|

**E\) Pair in reverse order, no overlap (can occur when the template fragment is a circular RNA)**


                                                                      start 1                         end 1


                                                                      -----------------read 1-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


       start 2                        end 2

       <---------------read 2--------------

       -----------------------------calculated fragment--------------------------------------->

Fragment building for different layouts of pairs that map to the reverse strand
-------------------------------------------------------------------------------

**A\) Pair in order, no overlap**


       start 2                         end 2

       -----------------read 2-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                                                      start 1                     end 1

                                                                      <---------------read 1------------

       <---------------------------calculated fragment---------------------------------------

|

|

|


**B\) Pair in order, overlap**


       start 2                         end 2

       -----------------read 2-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                            start 1                       end 1

                                            <---------------read 1---------------

       <-------------------calculated fragment------------------------

|

|

|

**C\) Pair in order, identical overlap**


                                            start 2                         end 2

                                            -----------------read 2---------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                            start 1                         end 1

                                            <---------------read 1-----------------

                                            <----calculated fragment-----

|

|

|

**D\) Pair in reverse order, both reads exceed each other (template fragment length is smaller than read length)**


                                                   start 2                         end 2

                                                   -----------------read 2-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


                                     start 1                         end 1

                                       <------------read 1---------------

    calculated fragment:              <------------------------

|

|

|


**E\) Pair in reverse order, no overlap (can occur when the template fragment is a circular RNA)**


                                                                      start 2                         end 2


                                                                      -----------------read 2-------------->


------------------reference sequence: forward strand------------------------------------------------------------------------------->


       start 1                         end 1

       <---------------read 1--------------

       <----------------------------calculated fragment----------------------------------------

