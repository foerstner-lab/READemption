Installation
============

Requirements
------------

TRAPL was developed using Python 3.3 and for best performance the user
is advised to run TRAPL with this version, too. Any other Python 3
version should work as well. Also Python 2.7 can be used if the
library `futures <https://pypi.python.org/pypi/futures>`_ is
installed. In any case, the third party modules `pysam
<https://code.google.com/p/pysam>`_ as well as `setuptool
<https://pypi.python.org/pypi/setuptools>`_ and `pip
<http://www.pip-installer.org>`_ in order to make the installation
easy by retrieving are required. TRAPL uses the short read mapper
`segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ for
the mapping and this software needs to be installed. The subcommand
`viz_align`, `viz_gene_quanti`, `viz_deseq` require the Python library
`Matplotlib <http://matplotlib.org/>`_. `R
<http://www.r-project.org/>`_ and the bioconductor package `DESeq
<http://bioconductor.org/packages/release/bioc/html/DESeq.html>`_ are
necessary for the subcommand `deseq` which performs differential gene
expression analysis. Don't worry - in the following the installation
of all these requirements will be covered.

Installing on a fresh Ubuntu image
----------------------------------

The following installation procedure was tested on a 
`Amazon AWS t1.micro
<http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts_micro_instances.html>`_
instance with Ubuntu Server 13.10 image.

Before starting it is a good idea to update the package list::

  sudo apt-get update

Ubuntu 13.10 has Python 3.3 already installed. If this is not the case
install::

 sudo apt-get install python3

Install setuptools::

 sudo apt-get install python3-setuptools

Install Matplotlib::

 sudo apt-get install python3-matplotlib

Additionally, Cython is needed.::

  sudo apt-get install cython3
  sudo apt-get install zlib1g-dev

If PIP is not yet install you should get this, too.::

  curl https://raw.github.com/pypa/pip/master/contrib/get-pip.py > get-pip.py
  sudo python3.3 get-pip.py

Now you can use PIP to install pysam and TRAPL::

  pip-3.3 install pysam
  pip-3.3 install trapl

Install make and ncurses dev library.::

  sudo apt-get install make
  sudo apt-get install libncurses5-dev

Install segemehl.::

  curl http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_1_6.tar.gz > segemehl_0_1_6.tar.gz
  tar xzf segemehl_0_1_6.tar.gz
  cd segemehl_*/segemehl/ && make && cd ../../

Copying it.::

  sudo cp segemehl*/segemehl/segemehl.x /usr/bin/segemehl

Alternative.::

  mkdir ~/bin
  cp segemehl_0_1_7/segemehl/segemehl.x ~/bin

Install R::

  sudo apt-get install r-base

and libxml2 which is required for the installation of some R-packages.::

 sudo apt-get install libxml2-dev

Install DESeq in ::

  echo 'source("http://bioconductor.org/biocLite.R");biocLite("DESeq")' | Rscript -


..
.. Global installation
.. -------------------
.. 
.. Installation in the home directory of the user
.. ----------------------------------------------
.. 
.. Installation in a pyvenv
.. ----------------------
