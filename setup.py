#!/usr/bin/env python3
"""Create setup script for instalation of embed fast5"""
########################################################################
# File: setup.py
#  executablgite: setup.py
#
# Author: Andrew Bailey
# History: 12/09/17 Created
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
        description='Python utility functions',
        url='https://github.com/adbailey4/python_utils',
        author='Andrew Bailey',
        license='MIT',
        author_email='andbaile@ucsc.com',
        tests_require=["pytest>=4.0.0"],
        install_requires=['py3helpers>=0.2.1',
                          'pandas=>=0.24.2',
                          'h5py=2.9.0'],
        zip_safe=True
    )

    stop = timer()
    print("Running Time = {} seconds".format(stop-start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
