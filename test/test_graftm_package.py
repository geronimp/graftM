#!/usr/bin/env python

#=======================================================================
# Authors: Ben Woodcroft, Joel Boyd
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import sys
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.graftm_package import GraftMPackage

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    def test_acquire(self):
        pkg = GraftMPackage.acquire(os.path.join(path_to_data, 'mcrA.gpkg'))
        self.assertEqual(os.path.join(path_to_data, 'mcrA.gpkg','mcrA.hmm'),
                         pkg.alignment_hmm_path())
        self.assertEqual(False, pkg.use_hmm_trusted_cutoff())
        self.assertEqual(os.path.join(path_to_data, 'mcrA.gpkg','mcrA.refpkg','treeoN87dL.tre'),
                         pkg.reference_package_tree_path())
        
    def test_version3_unaligned_fasta(self):
        pkg = GraftMPackage.acquire(os.path.join(path_to_data, 'mcrA.gpkg'))
        self.assertEqual(3, pkg.version)
        self.assertEqual(os.path.join(path_to_data, 'mcrA.gpkg','mcrA.faa'),
                         pkg.unaligned_sequence_database_path())

    def test_is_protein_package(self):
        pkg = GraftMPackage.acquire(os.path.join(path_to_data, 'mcrA.gpkg'))
        self.assertEqual(True, pkg.is_protein_package())
        self.assertEqual(True, pkg.is_protein_package())
        pkg = GraftMPackage.acquire(os.path.join(path_to_data, '61_otus.gpkg'))
        self.assertEqual(False, pkg.is_protein_package())
        self.assertEqual(False, pkg.is_protein_package())
        

if __name__ == "__main__":
    unittest.main()
