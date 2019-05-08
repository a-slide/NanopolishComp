# -*- coding: utf-8 -*-

# Define self package variable
__version__ = "0.6.2"
__description__ = 'NanopolishComp is a Python3 package for downstream analyses of Nanopolish output files'

# Collect info in a dictionary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __description__,
    "url": "https://github.com/a-slide/NanopolishComp",
    "author": 'Adrien Leger',
    "author_email": 'aleg@ebi.ac.uk',
    "license": "MIT",
    "python_requires":'>=3.3',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'],
    "install_requires": ['numpy>=1.14.0', 'tqdm>=4.23.4'],
    "packages": [__name__],
    "entry_points":{'console_scripts': [
        'NanopolishComp=NanopolishComp.__main__:main']}}
