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
import os
import sys
from string import split

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.deduplicator import Deduplicator
from graftm.sequence_io import Sequence

class Tests(unittest.TestCase):
    def s(self, seqs):
        to_return = []
        current_seq = None
        for bit in split(seqs,' '):
            if current_seq:
                current_seq.seq = bit
                to_return.append(current_seq)
                current_seq = None
            else:
                current_seq = Sequence(bit, '')
        return to_return
    
    def nd(self, seqs):
        return sorted([' '.join([x.name for x in s])
                       for s in self.deduplicator.deduplicate(seqs)])
        
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.deduplicator = Deduplicator()
        self.d = self.deduplicator
        
    def test_deduplicate_hello_world(self):
        seqs = self.s('1 AAA 2 AAT')
        self.assertEqual(['1','2'], self.nd(seqs))
                
    def test_deduplicate_one_dup(self):
        self.assertEqual(['1 2'], self.nd(self.s('1 AAA 2 AAA')))
                
    def test_deduplicate_one_dup_and_one_not(self):
        self.assertEqual(['1 3','2'], self.nd(self.s('1 AAA 2 ATA 3 AAA')))
        
    def test_lca(self):
        self.assertEqual(['1 3','2'], self.nd(self.s('1 AAA 2 ATA 3 AAA')))
        dees = self.d.deduplicate(self.s('1 AAA 2 ATA 3 AAA'))
        self.assertEqual([['one','two'], ['one','two','five']],
                         self.d.lca_taxonomy(dees, 
                                             {'1': ['one','two','four'],
                                              '2': ['one','two','five'],
                                              '3': ['one','two','three']}))
        
    def test_lca_more_specific_then_less(self):
        self.assertEqual(['1 3','2'], self.nd(self.s('1 AAA 2 ATA 3 AAA')))
        dees = self.d.deduplicate(self.s('1 AAA 2 ATA 3 AAA'))
        self.assertEqual([['one'], ['one','two']],
                         self.d.lca_taxonomy(dees, 
                                             {'1': ['one','two','four'],
                                              '2': ['one','two'],
                                              '3': ['one']}))
        
    def test_group_of_three_deduplication(self):
        dees = self.d.deduplicate(self.s('1 AAA 2 AAA 3 AAA'))
        self.assertEqual([['one','two','three']],
                         self.d.lca_taxonomy(dees, 
                                             {'1': ['one','two','three','1'],
                                              '2': ['one','two','three','2'],
                                              '3': ['one','two','three','3']}))
        
if __name__ == "__main__":
    unittest.main()
