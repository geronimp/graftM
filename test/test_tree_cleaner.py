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
from skbio.tree import TreeNode
from string import split as _

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.tree_cleaner import TreeCleaner

class Tests(unittest.TestCase):
    
    def match(self, tree, names):
        tc  = TreeCleaner()
        return tc.match_alignment_and_tree_sequence_ids(\
                                    names, TreeNode.read(StringIO(tree)))

    def test_match_alignment_and_tree_sequence_ids(self):
        self.match(u'(a,(b,c));',_('a b c'))
        
    def test_match_alignment_and_tree_sequence_ids_tree_not_align(self):
        self.assertRaises(Exception, self.match, u'(a,(b,c));',_('a b'))
        
    def test_match_alignment_and_tree_sequence_ids_align_not_tree(self):
        self.assertRaises(Exception, self.match, u'(a,(b,c));',_('a b c d'))
        
    def test_match_alignment_and_tree_sequence_ids_underscores(self):
        self.match(u'(\'a_2\',(b,c));',_('a_2 b c'))
        
    def test_remove_sequences_easy(self):
        tc  = TreeCleaner()
        tree = TreeNode.read(StringIO(u'(((a,b)F,c)d);'))
        tc.remove_sequences(tree, ['b'])
        self.assertEqual('((c,a)d);\n', str(tree))
        
    def test_remove_sequences_underscores(self):
        tc  = TreeCleaner()
        tree = TreeNode.read(StringIO(u"(((a,b),(c_yeh,L))d);"))
        tc.remove_sequences(tree, ['c_yeh'])
        self.assertEqual('(((a,b),L)d);\n', str(tree))
        
        
    def test_remove_sequences_with_named_internal_nodes(self):
        tc  = TreeCleaner()
        tree = TreeNode.read(StringIO(u"('Asulf_Archaeoglobus.1_2280~2522125074':7.17,(('Afulgi_764~2528311132':0.0,'CP006577_764~2588253768':0.0):0.0,'AE000782_746~638154502':0.0)'s__Archaeoglobus fulgidus':7.555):1.461;\n"))
        tc.remove_sequences(tree, 
                            ['CP006577_764~2588253768',
                             'Afulgi_764~2528311132'])
        self.assertEqual("('Asulf_Archaeoglobus.1_2280~2522125074':7.17,'AE000782_746~638154502':7.555):1.461;\n", str(tree))

if __name__ == "__main__":
    unittest.main()
