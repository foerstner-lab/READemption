try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='READemption',
    version='0.1.9',
    packages=['reademptionlib', 'tests'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='READemption - A RNA-Seq Analysis Pipeline',
    url='',
    install_requires=[
        "pysam >= 0.7.6"
    ],
    scripts=['bin/reademption'],
    license='LICENSE.txt',
    long_description=open('README.txt').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)'
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
