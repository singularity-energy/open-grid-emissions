#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="open-grid-emissions",
    # setuptools_scm lets us automagically get package version from GitHub tags
    setup_requires=["setuptools_scm"],
    use_scm_version=True,
    # Directory to search recursively for __init__.py files defining Python packages
    packages=find_packages("src"),
    # Location of the "root" package:
    package_dir={"": "src"},
)
