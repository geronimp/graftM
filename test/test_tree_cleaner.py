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
import tempfile
import sys

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','graftM')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.tree_cleaner import TreeCleaner

class Tests(unittest.TestCase):
    
    def match(self, tree, alignment):
        tc  = TreeCleaner()
        with tempfile.NamedTemporaryFile() as aln:
            aln.write(alignment)
            aln.flush()
            with tempfile.NamedTemporaryFile(suffix='.nwk') as tref:
                tref.write(tree)
                tref.flush()
                
                return tc.match_alignment_and_tree_sequence_ids(aln.name,
                                                              tref.name)

    def test_match_alignment_and_tree_sequence_ids(self):
        self.match('(a,(b,c));',"\n".join('>a A >b A >c C'.split(' ')))
        
    def test_match_alignment_and_tree_sequence_ids_tree_not_align(self):
        self.assertRaises(Exception, self.match, '(a,(b,c));',"\n".join('>a A >b A'.split(' ')))
        
    def test_match_alignment_and_tree_sequence_ids_align_not_tree(self):
        self.assertRaises(Exception, self.match, '(a,(b,c));',"\n".join('>a A >b A >c C >d T'.split(' ')))
        
    def test_match_alignment_and_tree_sequence_ids_underscores(self):
        self.match('(a_2,(b,c));',"\n".join('>a_2 A >b A >c C'.split(' ')))
        
    def test_match_alignment_and_tree_sequence_ids_comments(self):
        self.match('(a,(b,c));',"\n".join('>a comment A >b A >c C'.split(' ')))
        

if __name__ == "__main__":
    unittest.main()
