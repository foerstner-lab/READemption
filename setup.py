try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='READemption',
    version='0.3.5',
    packages=['reademptionlib', 'tests'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='A RNA-Seq Analysis Pipeline',
    url='',
    install_requires=[
        "biopython >= 1.65",
        "matplotlib >= 1.4.2",
        "pandas >= 0.15.2",
        "pysam >= 0.8.1"
    ],
    scripts=['bin/reademption'],
    license='ISC License (ISCL)',
    long_description=open('README.rst').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
