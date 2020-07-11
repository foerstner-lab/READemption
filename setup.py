try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="READemption",
    version="1.0.1",
    packages=["reademptionlib", "tests"],
    author="Konrad U. Förstner, Till Sauerwein",
    author_email="konrad@foerstner.org",
    description="A RNA-Seq Analysis Pipeline",
    url="",
    install_requires=[
        "biopython >= 1.73",
        "matplotlib >= 2.2.2",
        "pandas >= 0.24.2",
        "pysam >= 0.15.2",
        "sphinx-argparse >=0.2.5"
    ],
    scripts=["bin/reademption"],
    license="ISC License (ISCL)",
    long_description=open("README.rst").read(),
    classifiers=[
        "License :: OSI Approved :: ISC License (ISCL)",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
