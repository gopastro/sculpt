#!/usr/bin/env python
from setuptools import setup, find_packages

NAME = 'sculpt'
VERSION = '0.1.dev'

setup(
    name=NAME,
    version=VERSION,
    description='Spectral Line Data Cube Reduction and Visualization Tool for Radio Astronomy',
    author='Sculpt Developers',
    packages=find_packages(),
    )
