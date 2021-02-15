#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on

@author: Lilian, Jorge, Aitor
"""

from setuptools import setup
from os import path

#this_directory = path.abspath(path.dirname(__file__))
#with open(path.join(this_directory, 'theory.md'), encoding='utf-8') as f: this is the theory file that need the be delivered.
#    long_description_in = f.read()

setup(name='promod',
      version='0',
      description='A macromolecular complex builder',
      long_description="",
      long_description_content_type='text/markdown',
      author='Jorge Dominguez, Lilian Marie Boll, Aitor Gonzalez',
      author_email='lilianmarie.boll01@estudiant.upf.edu, aitor.gonzalez03@estudiant.upf.edu, jorge.dominguez01@estudiant.upf.edu',
      license='MIT',
      keywords='bioinformatics structural_bioinformatics protein complex builder modeler',
      packages=['builder'],
      scripts=['promod', 'promod-tk','scripts/pairpdbs.py', 'scripts/pdbsplit.py'],
      install_requieres = ['Biopython'],
      url='https://github.com/jorgedb98/SBI',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3',
      include_package_data=True)
