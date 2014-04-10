Troubleshooting
===============

Pickling error
--------------

Problem
~~~~~~~

When running reademption (e.g. with Python 3.2) you get an error like this:
::
   Traceback (most recent call last):
     File "/usr/lib/python3.2/multiprocessing/queues.py", line 272, in _feed
      send(obj)
   _pickle.PicklingError: Can't pickle <class 'method'>: attribute lookup builtins.method failed

Solution
~~~~~~~~

Switch to Python 3.3 or higher. The newer versions can handle this.

DESeq2 unused argument error
----------------------------

Problem
~~~~~~~

Whenn running ``deseq`` you get an error similar to this:
::
   Error in results(dds, contrast = c("condition", "XXX", "YYY")) : 
     unused argument (contrast = c("condition", "XXX", "YYY"))
   Execution halted
   Apparently DESeq did not generate the file "deseq_comp_XXX_vs_YYY.csv". Extension stopped.

Solution
~~~~~~~~

Get a newer version of DESeq2. Version 1.2.10 and higher should work.
