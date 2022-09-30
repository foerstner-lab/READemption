try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="READemption",
    version="2.0.3",
    packages=["reademptionlib", "tests"],
    author="Konrad U. Förstner, Till Sauerwein",
    author_email="konrad@foerstner.org",
    description="A RNA-Seq Analysis Pipeline",
    url="",
    install_requires=[
        "biopython >= 1.79",
        "matplotlib >= 3.5.2",
        "pandas >= 1.4.3",
        "pysam >= 0.19.1",
        "seaborn >= 0.11.2",
        "sphinx-argparse >=0.2.5"
    ],
    scripts=["bin/reademption"],
    license="MIT License",
    long_description=open("README.rst").read(),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
