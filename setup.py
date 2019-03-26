#!/usr/bin/env python3
"""Create setup script for instalation of embed fast5"""
########################################################################
# File: setup.py
#  executablgite: setup.py
#
# Author: Andrew Bailey
# History: 3/21/19 Created
########################################################################

import sys
from timeit import default_timer as timer
from setuptools import setup, find_packages


def main():
    """Main docstring"""
    start = timer()
    setup(
        name="embed",
        version='0.0.3',
        description='Embed fast5 with event table from nanopolish model',
        url='https://github.com/adbailey4/embed_fast5',
        author='Andrew Bailey',
        license='MIT',
        packages=['embed'],
        author_email='andbaile@ucsc.com',
        install_requires=['py3helpers>=0.2.1',
                          'pandas>=0.24.2',
                          'h5py>=2.9.0'],
        zip_safe=True
    )

    stop = timer()
    print("Running Time = {} seconds".format(stop-start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
