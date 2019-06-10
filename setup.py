#!/usr/bin/env python3
"""Create setup script for instalation of embed fast5"""
########################################################################
# File: setup.py
#  executablgite: setup.py
#
# Author: Andrew Bailey
# History: 3/21/19 Created
########################################################################
from timeit import default_timer as timer

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j8']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='embed',
    version='0.0.1',
    author='Dean Moldovan',
    author_email='dean0x7d@gmail.com',
    description='A test project using pybind11 and CMake',
    long_description='',
    ext_modules=[CMakeExtension('cmake_example3')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)


# def main():
#     """Main docstring"""
#     start = timer()
#     setup(
#         name="embed",
#         version='0.0.5',
#         description='Embed fast5 with event table from nanopolish model',
#         url='https://github.com/adbailey4/embed_fast5',
#         author='Andrew Bailey',
#         license='MIT',
#         packages=['embed'],
#         ext_modules=[CMakeExtension('cmake_example2')],
#         cmdclass=dict(build_ext=CMakeBuild),
#         # scripts=["bin/split_multi_fast5", "bin/embed_fast5s"],
#         author_email='andbaile@ucsc.com',
#         # install_requires=['py3helpers[seq_tools]>=0.2.9',
#         #                   'pandas>=0.24.2',
#         #                   'h5py>=2.9.0'],
#         zip_safe=False
#     )
#
#     stop = timer()
#     print("Running Time = {} seconds".format(stop-start), file=sys.stderr)
#
#
# if __name__ == "__main__":
#     main()
#     raise SystemExit
