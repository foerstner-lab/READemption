About
-----

READemption is a pipeline for the computational evaluation of RNA-Seq
data. It was originally developed to process dRNA-Seq reads (as
introduced by Sharma et al., Nature, 2010) originating from bacterial
samples. Meanwhile is has been extended to process data generated in
different experimental setups and from all domains of life. The
functions which are accessible via a command-line interface cover read
processing and aligning, coverage calculation, gene expression
quantification, differential gene expression analysis as well as
visualization. In order to set up and perform analyses quickly
READemption follows the principal of "convention over configuration":
Once the input files are copied/linked into defined folders no further
parameters have to be given. Still, READemption's behavior can be
adapted to specific needs of the user by parameters.

Documentation
-------------

Documentation can be found on
`here <http://pythonhosted.org/READemption/>`__.

Installation
------------

Short version (if you have all the requirements installed):

::

    $ pip install READemption

`Long version <http://pythonhosted.org/READemption/installation.html>`__
including a decryption of the requirements and how do you get them.

License
-------

`ICSL <https://en.wikipedia.org/wiki/ISC_license>`__ (Internet Systems
Consortium license ~ simplified BSD license) - see LICENSE.txt

Development
-----------

-  If possible follow the principal of "convention over configuration".
   This means input file are copied/linked into a fixed location and the
   resulting files are placed in fixed locations.

-  The classes should be path agnostic as far a possible. The controller
   is taking care of that and calls them adequately.

-  The git braching model is very close to the one proposed
   `here <http://nvie.com/posts/a-successful-git-branching-model/>`__.
   There two main branches:

   -  master
   -  dev(elopment)

   And there are further supporting branches:

   -  feature branches - branched off and back to the dev branch
   -  release branches - branched off from dev and merged back into dev
      and master
   -  hotfix branches - branched off from master and merged back into
      dev and master


