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
import tempfile
import os
import sys
import io
from biom.util import biom_open
import subprocess

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.summarise import Stats_And_Summary

class Tests(unittest.TestCase):
    def test_iterate_otu_table_rows_hello_world(self):
        s = Stats_And_Summary()
        self.assertEqual(
                         [(1, ['ab','c'], [1])],
                         list(s._iterate_otu_table_rows([{'readname': ['ab','c']}]))
                         )  
        
    def test_iterate_otu_table_rows_two_samples_same_tax(self):
        s = Stats_And_Summary()
        self.assertEqual(
                         [(1, ['ab','c'], [1,1])],
                         list(s._iterate_otu_table_rows(({'readname': ['ab','c']},
                                                         {'readname2': ['ab','c']})))
                         )  
        
    def test_iterate_otu_table_rows_two_samples_different_counts(self):
        s = Stats_And_Summary()
        self.assertEqual(
                         [(1, ['ab','c'], [2,1])],
                         list(s._iterate_otu_table_rows([
                                                         {'readname': ['ab','c'], 'readname23': ['ab','c']},
                                                         {'readname2': ['ab','c']}
                                                         ]))
                         )  
        
    def test_iterate_otu_table_rows_two_samples_different_counts_twotax(self):
        s = Stats_And_Summary()
        self.assertEqual(
                         [(1, ['ab','d'], [1,0]), (2, ['ab','c'], [1,1])],
                         list(s._iterate_otu_table_rows([
                                                         {'readname': ['ab','c'], 'readname23': ['ab','d']},
                                                         {'readname2': ['ab','c']}
                                                         ]))
                         )
        
    def test_iterate_otu_table_rows_two_samples_second_new_tax(self):
        s = Stats_And_Summary()
        self.assertEqual(
                         [(1, ['ab','e'], [0,1]), (2, ['ab','c'], [2,0])],
                         list(s._iterate_otu_table_rows([
                                                         {'readname': ['ab','c'], 'readname23': ['ab','c']},
                                                         {'readname2': ['ab','e']}
                                                         ]))
                         )
        
    def test_write_otu_table(self):
        string = io.StringIO()
        s = Stats_And_Summary()
        s.write_tabular_otu_table(('sample1','sample2'),
                                  [
                                                         {'readname': ['ab','c']},
                                                         {'readname2': ['ab','c']}
                                                         ],
                                  string)
        self.assertEqual(u'#ID\tsample1\tsample2\tConsensusLineage\n1\t1\t1\tab; c\n', string.getvalue())
        
    def test_write_biom(self):
        with tempfile.NamedTemporaryFile(suffix='biom') as biom:
            with biom_open(biom.name,'w') as f:
                s = Stats_And_Summary()
                s.write_biom(('sample1','sample2'),
                          [
                                                 {'readname': ['ab','c'], 'readnameE': ['ab','d']},
                                                 {'readname2': ['ab','c']}
                                                 ],
                          f)

            with tempfile.NamedTemporaryFile(suffix='csv') as biom_out:
                os.remove(biom_out.name) #delete because otherwise biom complains
                subprocess.check_call("biom convert -i %s -o %s --table-type 'OTU table' --to-tsv --header-key taxonomy" %
                                                   (biom.name, biom_out.name), shell=True)
                observed = open(biom_out.name).read()
                self.assertEqual('''# Constructed from biom file
#OTU ID\tsample1\tsample2\ttaxonomy
1\t1.0\t0.0\tab; d
2\t1.0\t1.0\tab; c
''', observed)
        
        

if __name__ == "__main__":
    unittest.main()
