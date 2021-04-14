#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 2021

@author: Lilian, Jorge, Aitor
"""

from setuptools import setup
from os import path

with open('doc.md','r') as f:
    long_description=f.read()

setup(name='complex_asssembler',
      version='0.1',
      description='A macromolecular complex builder',
      long_description_content_type='text/markdown',
      author='Jorge Dominguez, Lilian Marie Boll, Aitor Gonzalez',
      author_email='jorge.dominguez01@estudiant.upf.edu, lilianmarie.boll01@estudiant.upf.edu, aitor.gonzalez03@estudiant.upf.edu',
      license='MIT',
      keywords='bioinformatics structural_bioinformatics protein complex builder modeler',
      packages=['complex_asssembler'],
      scripts=['complex_assembler.py', 'functions.py'],
      install_requires = ['Biopython'],
      url='https://github.com/jorgedb98/SBI',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3',
      include_package_data=True)
