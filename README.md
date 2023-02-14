[![Latest Version](https://img.shields.io/pypi/v/reademption.svg)](https://pypi.python.org/pypi/READemption/)
[![License](https://img.shields.io/pypi/l/reademption.svg)](https://pypi.python.org/pypi/READemption/)
![Build Status](https://github.com/foerstner-lab/READemption/actions/workflows/main.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/reademption/badge/?version=latest)](https://reademption.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/18210971.svg)](https://zenodo.org/badge/latestdoi/18210971)
[![Anaconda-Server Badge](https://anaconda.org/till_sauerwein/reademption/badges/version.svg)](https://conda.anaconda.org/till_sauerwein)
<img src="./coverage.svg">

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

Documentation can be found on [here](https://reademption.readthedocs.io).

Installation
------------

Short version (if you have all the requirements installed):

    $ pip install READemption

[Long version](https://reademption.readthedocs.io)
including a description of the requirements and how do you get them.

License
-------

[MIT](https://opensource.org/licenses/MIT) 
(MIT License) - see LICENSE.txt

Development
-----------

* If possible follow the principal of "convention over
  configuration". This means input file are copied/linked into a fixed
  location and the resulting files are placed in fixed locations.

* The classes should be path agnostic as far a possible. The controller
  is taking care of that and calls them adequately.

* The git braching model is very close to the one 
  proposed [here](http://nvie.com/posts/a-successful-git-branching-model/).
  There two main branches:
    * master 
    * dev(elopment)

    And there are further supporting branches:
    * feature branches - branched off and back to the dev branch
    * release branches - branched off from dev and merged back into
                       dev and master
    * hotfix branches - branched off from master and merged back into
                      dev and master
