#!/usr/bin/env python3
"""Create setup script for installation of embed_fast5"""
########################################################################
# File: setup.py
#  executable: setup.py
#
# Author: Andrew Bailey (https://www.benjack.io/2018/02/02/python-cpp-revisited.html)
# History: 3/21/19 Created
########################################################################

import os
import re
import sys
import platform
import subprocess
from shutil import copyfile, copymode
from timeit import default_timer as timer

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


def copy_test_file(src_file):
    """
    Copy ``src_file`` to `tests/bin` directory, ensuring parent directory
    exists. Messages like `creating directory /path/to/package` and
    `copying directory /src/path/to/package -> path/to/package` are
    displayed on standard output. Adapted from scikit-build.
    """
    # Create directory if needed
    dest_dir = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), 'tests', 'bin')
    if dest_dir != "" and not os.path.exists(dest_dir):
        print("creating directory {}".format(dest_dir))
        os.makedirs(dest_dir)

    # Copy file
    dest_file = os.path.join(dest_dir, os.path.basename(src_file))
    print("copying {} -> {}".format(src_file, dest_file))
    copyfile(src_file, dest_file)
    copymode(src_file, dest_file)


def copy_test_file2(src_file, dest_dir):
    """
    Copy ``src_file`` to `dest_dir` directory, ensuring parent directory
    exists. Messages like `creating directory /path/to/package` and
    `copying directory /src/path/to/package -> path/to/package` are
    displayed on standard output. Adapted from scikit-build.
    """
    # Create directory if needed
    if dest_dir != "" and not os.path.exists(dest_dir):
        print("creating directory {}".format(dest_dir))
        os.makedirs(dest_dir)

    # Copy file
    dest_file = os.path.join(dest_dir, os.path.basename(src_file))
    print("copying {} -> {}".format(src_file, dest_file))
    copyfile(src_file, dest_file)
    copymode(src_file, dest_file)


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
        shared_opt = "OFF"
        if "BUILD_SHARED_LIBS" in os.environ:
            shared_opt = os.environ["BUILD_SHARED_LIBS"]
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DBUILD_SHARED_LIBS={}'.format(shared_opt),
                      '-DCMAKE_VERBOSE_MAKEFILE=ON']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            # build_args += ['--', '-j8']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        # subprocess.check_call(['cmake', '--install', '.'], cwd=self.build_temp)
        # Copy *_test file to tests directory
        embed_cpp_tests = os.path.join(self.build_temp, "tests", 'test_embed')
        embed_main = os.path.join(self.build_temp, 'embed_main')
        copy_test_file(embed_cpp_tests)
        copy_test_file(embed_main)
        print()  # Add empty line for nicer output


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        install.run(self)
        build_temp = self.build_lib.replace("lib", "temp")
        source = os.path.join(os.path.dirname(os.path.abspath(__file__)), build_temp, "embed_main")
        target = os.path.join(self.install_scripts, "embed_main")
        if os.path.isfile(target):
            os.remove(target)
        self.copy_file(source, target)

        source = os.path.join(os.path.dirname(os.path.abspath(__file__)), build_temp, "tests/test_embed")
        target = os.path.join(self.install_scripts, "test_embed")
        if os.path.isfile(target):
            os.remove(target)

        self.copy_file(source, target)


def get_version():
    try:
        content = open("CMakeLists.txt", "r").read()
        version = re.search(r'project\(embed_fast5 VERSION (.*)\)', content).group(1)
        return version.strip()
    except RuntimeError:
        return None


def main():
    """Main docstring"""
    start = timer()
    setup(
        name="embed",
        version=get_version(),
        description='Embed fast5 with event table from nanopolish model',
        url="https://github.com/adbailey4/embed_fast5",
        author='Andrew Bailey',
        license='MIT',
        packages=find_packages('src'),
        package_dir={'': 'src'},
        ext_modules=[CMakeExtension('embed.bindings')],
        cmdclass=dict(build_ext=CMakeBuild, install=PostInstallCommand),
        scripts=["src/scripts/split_multi_fast5.py", "src/scripts/embed_fast5s.py"],
        author_email='andbaile@ucsc.com',
        install_requires=['py3helpers>=1.0.0',
                          'h5py>=2.9.0'],
        zip_safe=False,
        test_suite='tests'
    )

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
