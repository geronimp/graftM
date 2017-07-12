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
import os.path
import sys
from StringIO import StringIO
import tempfile

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.greengenes_taxonomy import GreenGenesTaxonomy,\
    MalformedGreenGenesTaxonomyException, DuplicateTaxonomyException

class Tests(unittest.TestCase):
    def test_read_hello_world(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria']},\
                          GreenGenesTaxonomy.read(StringIO('seq1\tbacteria; cyanobacteria')).taxonomy)

    def test_read_semicolon_no_space(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria']},\
                          GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria')).taxonomy)

    def test_raises_when_incorrect_num_fields(self):
        with self.assertRaises(MalformedGreenGenesTaxonomyException):
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria\n'\
                                           'seq2\n'
                                           ))

    def test_ok_when_taxonomy_empty(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria'], 'seq2': []},\
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria\n'\
                                           'seq2\t\n'
                                           )).taxonomy)

    def test_raises_when_duplicate_names(self):
        with self.assertRaises(DuplicateTaxonomyException):
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria\n'\
                                           'seq1\tbacteria;cyanobacteria\n'
                                           ))

    def test_raises_when_missing_middle(self):
        with self.assertRaises(MalformedGreenGenesTaxonomyException):
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria\n'\
                                           'seq2\tbacteria;;cyanobacteria\n'
                                           ))

    def test_removes_empties_at_end(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria'], 'seq2': ['bacteria','bluebacteria']},\
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria;\n'\
                                           'seq2\tbacteria;bluebacteria;;\n'
                                           )).taxonomy)

    def test_ignores_empty_lines(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria'], 'seq2': ['bacteria','bluebacteria']},\
            GreenGenesTaxonomy.read(StringIO('seq1\tbacteria;cyanobacteria;\n'\
                                           'seq2\tbacteria;bluebacteria;;\n'\
                                           '\n'
                                           )).taxonomy)

    def test_strip_identifier(self):
        self.assertEqual({'seq1': ['bacteria','cyanobacteria'], 'seq2': ['bacteria','bluebacteria']},\
            GreenGenesTaxonomy.read(StringIO('seq1 \tbacteria;cyanobacteria;\n'\
                                           'seq2\tbacteria;bluebacteria;;\n'
                                           )).taxonomy)

    def test_input_file(self):
        with tempfile.NamedTemporaryFile(prefix='graftm_greengenes_tax_testing') as tf:
            tf.write('seq1\tbacteria;cyanobacteria')
            tf.flush()
            self.assertEqual({'seq1': ['bacteria','cyanobacteria']},\
                      GreenGenesTaxonomy.read_file(tf.name).taxonomy)

if __name__ == "__main__":
    unittest.main()
