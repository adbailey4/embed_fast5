#!/usr/bin/env python3
"""Testing script wrapper for cpp tests executed by 'python setup.py test' """
########################################################################
# File: test_cpp_tests.py
#  executable: test_cpp_tests.py
#
# Author: Andrew Bailey (https://www.benjack.io/2018/02/02/python-cpp-revisited.html)
# History: 08/21/19 Created
########################################################################


import unittest
import subprocess
import os


class CppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(CppTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])

    def test_cpp(self):
        print("\n\nTesting C++ code...")
        test_path = os.path.join(self.HOME, 'tests/bin', 'test_embed')
        print()
        subprocess.check_call([test_path, self.HOME])
        print()
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
