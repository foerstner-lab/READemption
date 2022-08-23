Installation and updating
=========================

Requirements
------------

READemption was started to be developed using Python 3.2 and the user
is advised to run READemption with this or a higher version. In any
case `setuptools <https://pypi.python.org/pypi/setuptools>`_ and `pip
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

Installing on a fresh Ubuntu system
-----------------------------------


1. Installing all required Debian/Ubuntu packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting it is a good idea to update the package list::

  sudo apt-get update

Now you can install the packages::

  sudo apt-get install python3 python3-setuptools python3-pip python3-matplotlib cython3 zlib1g-dev  make libncurses5-dev r-base libxml2-dev

Some comments:

- Python version 3.6 and higher is recommended
- ``cython`` is required for ``pysam``
- ``make``, ``libncurses5-dev`` and ``zlib1g-dev`` are needed for ``segemehl``
- ``libxml2`` is required for the installation of some of the R-packages

2. Install segemehl
~~~~~~~~~~~~~~~~~~~

::

  curl https://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/segemehl-0.3.4.tar.gz > segemehl-0.3.4.tar.gz
  tar segemehl-0.3.4.tar.gz
  cd segemehl_*/segemehl*/ && make all && cd ../../


Copying the executable to a location that is part of the ``PATH`` e.g
``/usr/bin/`` ...

::
  sudo cp segemehl-0.3.4/segemehl-0.3.4/segemehl.x /usr/bin/segemehl.x

... or the bin folder of your home directory::

  mkdir ~/bin
  cp segemehl-0.3.4/segemehl-0.3.4/segemehl.x ~/bin


Alternatively, segemehl can be installed via conda::

  conda install -c bioconda segemehl


3. Install DESeq2
~~~~~~~~~~~~~~~~~

Start ``R``::

  R


and install the DESeq2 package inside of the interactive command line
interface. You might be asked to confirm the installation path::

  source("http://bioconductor.org/biocLite.R")
  biocLite("DESeq2")

Leave ``R``::

  quit(save = "no")


Install  READemption and it's dependcies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now you can use ``pip`` to install ``READemption`` and further python
packages::

  sudo pip3 install READemption

Voilà! You should now be able to call READemption::

  reademption -h


Installing on a Apple OS X
--------------------------

(Many thanks to Lei Li for contribution this part!)

1. Installing all required software/packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To download and install Python 3 follow the instruction at this
`download page <https://www.python.org/downloads/>`_.

Download and install `xcode <https://developer.apple.com/xcode/>`_ (`page <https://developer.apple.com/xcode/downloads/>`_) and R
(download links are on `the frontpage <http://www.r-project.org/>`_).

To install ``pip`` open a terminal and run

::

  curl -O https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py python3 ez_setup.py # download and install pip 
  curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py 
  python3 get-pip.py

Install ``matplotlib``:

::

  pip3 install matplotlib


2. Installing segemehl, DESeq, pysam and READemption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The remaining installation steps are the same as descibed above. Just
open a terminal and run the commands.


Updating READemption
--------------------

Once you have installed READemption as described above you can easily
upgrade it to the newest version by running

::

  pip3 install --upgrade READemption


Package version requirements and installation
---------------------------------------------
Create a python version 3.9 environment
::
    conda create --name py3.9_reademption python=3.9

Activate the environment
::
    conda activate py3.9_reademption

Install pysam
::
    conda install -c bioconda pysam=0.19.1


Install segemehl
::
    conda install -c bioconda segemehl=0.3.4

Install biopython (numpy 1.23 included)

::
    conda install biopython=1.79

Install pandas

::
    conda install pandas=1.4.3

Install seaborn

::
    conda install seaborn=0.11.2

Install matplotlib

::
    conda install matplotlib=3.5.2

Install pytest (only required for testing)

::
    conda install pytest