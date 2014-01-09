try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='trapl',
    version='0.1.7',
    packages=['trapllib', 'tests'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='The RNA-Seq Analysis Pipeline',
    url='',
    install_requires=[
        "pysam >= 0.7.6"
    ],
    scripts=['bin/trapl'],
    license='LICENSE.txt',
    long_description=open('README.txt').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License',
        'Programming Language :: Python :: 3',
    ]
)
