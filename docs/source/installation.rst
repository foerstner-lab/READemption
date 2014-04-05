Installation
============

Requirements
------------

READemption was developed using Python 3.3 and the user is advised to
run READemption with this or a higher version. In any case, the third
party packages `pysam <https://code.google.com/p/pysam>`_ as well as
`setuptools <https://pypi.python.org/pypi/setuptools>`_ and `pip
<http://www.pip-installer.org>`_ should be available on the system in
order to make the installation easy. READemption uses the short read
mapper `segemehl
<http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ for the
mapping and this software needs to be installed. The subcommands
``viz_align``, ``viz_gene_quanti``, ``viz_deseq`` require the Python
library `Matplotlib <http://matplotlib.org/>`_. `R
<http://www.r-project.org/>`_ and the bioconductor package `DESeq2
<http://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ are
necessary for the subcommand ``deseq`` which performs differential
gene expression analysis. Don't worry - in the following the
installation of all these requirements will be covered.

Installing on a fresh Ubuntu installation
-----------------------------------------

The following installation procedure was tested on a `Amazon AWS
t1.micro
<http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts_micro_instances.html>`_
instance with Ubuntu Server 13.10 image.


1. Installing all required Debian/Ubuntu packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting it is a good idea to update the package list::

  sudo apt-get update

Now you can install the packages::

  sudo apt-get install python3 python3-setuptools python3-pip python3-matplotlib cython3 zlib1g-dev  make libncurses5-dev r-base libxml2-dev

Some comments:

- Ubuntu 13.10 should have Python 3.3 already installed.
- ``cython`` is required for ``pysam``
- ``make``, ``libncurses5-dev`` and ``zlib1g-dev`` are needed for ``segemehl``
- ``libxml2`` required for the installation of some of the R-packages

2. Install segemehl
~~~~~~~~~~~~~~~~~~~

::

  curl http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_1_7.tar.gz > segemehl_0_1_7.tar.gz
  tar xzf segemehl_0_1_7.tar.gz
  cd segemehl_*/segemehl/ && make && cd ../../

Copying the executable to a location that is part of the ``PATH`` e.g
``/usr/bin/`` ...

::

  sudo cp segemehl_0_1_7/segemehl/segemehl.x /usr/bin/segemehl.x
  sudo cp segemehl_0_1_7/segemehl/lack.x /usr/bin/lack.x

... or the bin folder of your home directory::

  mkdir ~/bin
  cp segemehl_0_1_7/segemehl/segemehl.x ~/bin

3. Install DESeq2
~~~~~~~~~~~~~~~~~

::

  echo 'source("http://bioconductor.org/biocLite.R");biocLite("DESeq2")' | sudo Rscript -

Install pysam and READemption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now you can use ``pip`` to install ``pysam`` and ``READemption``::

  sudo pip3 install pysam
  sudo pip3 install READemption

Voilà! You should now be able to call READemption::

  reademption -h

..
.. Global installation
.. -------------------
.. 
.. Installation in the home directory of the user
.. ----------------------------------------------
.. 
.. Installation in a pyvenv
.. ----------------------
