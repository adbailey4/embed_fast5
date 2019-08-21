#!/usr/bin/env python3
"""
    Place unit tests for python wrapped c++ functions
"""
########################################################################
# File: test_embedding_helpers.py
#  executable: test_embedding_helpers.py
#
# Author: Andrew Bailey
# History: 07/30/2019 Created
########################################################################

import unittest
import os
from embed import bindings


class TestPythonWrappers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestPythonWrappers, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-3])

    def test_add(self):
        self.assertEqual(3, bindings.add(1, 2))


if __name__ == '__main__':
    unittest.main()
