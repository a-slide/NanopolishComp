#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "NanopolishComp"
version = "0.6.4"
description = 'NanopolishComp is a Python3 package for downstream analyses of Nanopolish output files'
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/a-slide/NanopolishComp",
    author = 'Adrien Leger',
    author_email = 'aleg@ebi.ac.uk',
    license= "MIT",
    python_requires ='>=3.6',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'],
    install_requires = ['numpy>=1.14.0', 'tqdm>=4.23.4'],
    packages = [name],
    entry_points = {'console_scripts': ['NanopolishComp=NanopolishComp.__main__:main']})
