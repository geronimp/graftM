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
from io import StringIO
from dendropy import Tree

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.taxonomy_extractor import TaxonomyExtractor

class Tests(unittest.TestCase):
    
    def test_hello_world(self):
        self.assertEquals({u'a': [],
                           u'b': [u't1'],
                           u'c': [u't1']},
                          TaxonomyExtractor().taxonomy_from_annotated_tree(\
                            Tree.get(data="(a,(b,c)'t1':0.9)root;", schema='newick')))
        
    def test_bootstraps_in_annotated_tree(self):
        self.assertEquals({u'a': [],
                           u'b': [u't1'],
                           u'c': [u't1']},
                          TaxonomyExtractor().taxonomy_from_annotated_tree(\
                            Tree.get(data="(a,(b,c)'0.01973:t1':0.9)root;", schema='newick')))
        
        
    def test_bootstraps_in_annotated_tree_alongside_empty_taxa(self):
        self.assertEquals({u'a': [],
                           u'b': [],
                           u'c': ['tax'],
                           u'd': ['tax']},
                          TaxonomyExtractor().taxonomy_from_annotated_tree(\
                            Tree.get(data="(a,(b,(c,d:0.2)'0.2:tax')0.01973:0.9)root;", schema='newick')))
        
    def test_branch_lengths(self):
        '''https://github.com/geronimp/graftM/issues/192'''
        taxes = TaxonomyExtractor().taxonomy_from_annotated_tree(
            Tree.get(path=os.path.join(path_to_data, 'create', 'sulfitereductase.ben.tree'), schema='newick'))
        self.assertEquals([u'Aanerobic sulfite reductase asrC',
                           u'Anaerobic sulfite reductase asrC Group 3',
                           u'Unknown alpha and beta subunits',
                           u'0.856_PFAM_NIR_SIR,NIR_SIR_ferr'], # number is actually in the clade name
                          taxes['T506DRAFT_scaffold00010.10_60~2561511230'])
        

if __name__ == "__main__":
    unittest.main()
