#!/usr/bin/env python3

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
import tempfile
from io import StringIO
from dendropy import Tree

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.dendropy_tree_cleaner import DendropyTreeCleaner

class Tests(unittest.TestCase):
    def clean(self, cleaner, tree):
        with tempfile.NamedTemporaryFile() as f:
            cleaner.write_fasttree_newick(tree, f)
            f.flush()
            with open(f.name) as g:
                return g.read()

    def test_write_fasttree_newick(self):
        tc  = DendropyTreeCleaner()
        tree = Tree.get(data="((a,b),(d,e))root;", schema='newick')
        self.assertEqual("((a,b),(d,e));\n", self.clean(tc, tree))

        # Internal labels should be removed.
        tree = Tree.get(data="((a_2,b)c,(d,e)f)root;", schema='newick')
        self.assertEqual("((a_2,b),(d,e));\n", self.clean(tc, tree))

        # Quoted spaces should become underscores.
        tree = Tree.get(data="(('a 2',b),(d,e))root;", schema='newick')
        self.assertEqual("((a_2,b),(d,e));\n", self.clean(tc, tree))

        # Test underscores that are quoted.
        tree = Tree.get(data="(('a_2',b),(d,e))root;", schema='newick')
        self.assertEqual("((a_2,b),(d,e));\n", self.clean(tc, tree))

        # Test dashes
        tree = Tree.get(data="((ANME-2dV10_01644,b),(d,e))root;", schema='newick')
        self.assertEqual("((ANME-2dV10_01644,b),(d,e));\n", self.clean(tc, tree))

        # A more real world example with '~' characters (which never mattered actually).
        tree = Tree.get(
            data="('Asulf_Archaeoglobus.1_2280~2522125074':7.17,(('Afulgi_764~2528311132':0.0,'CP006577_764~2588253768':0.0):0.0,'AE000782_746~638154502':0.0)'s__Archaeoglobus fulgidus':7.555):1.461;\n",
            schema='newick')
        self.assertEqual("(Asulf_Archaeoglobus.1_2280~2522125074:7.17,((Afulgi_764~2528311132:0.0,CP006577_764~2588253768:0.0):0.0,AE000782_746~638154502:0.0):7.555):1.461;\n", self.clean(tc, tree))

    def match(self, tree, names):
        tc  = DendropyTreeCleaner()
        return tc.match_alignment_and_tree_sequence_ids(
            names, Tree.get(data=tree, schema='newick'))

    def test_match_alignment_and_tree_sequence_ids(self):
        self.match('(a,(b,c));',str.split('a b c'))

    def test_match_alignment_and_tree_sequence_ids_tree_not_align(self):
        self.assertRaises(Exception, self.match, '(a,(b,c));',str.split('a b'))

    def test_match_alignment_and_tree_sequence_ids_align_not_tree(self):
        self.assertRaises(Exception, self.match, '(a,(b,c));',str.split('a b c d'))

    def test_match_alignment_and_tree_sequence_ids_underscores(self):
        self.match('(\'a_2\',(b,c));',str.split('a_2 b c'))

    def test_remove_sequences_underscores(self):
        tc  = DendropyTreeCleaner()
        tree = Tree.get(data="(((a,b),(c_yeh,L))d);", schema='newick')
        tc.remove_sequences(tree, ['c_yeh'])
        self.assertEqual('((a,b),L)d', str(tree))

    def test_remove_sequences_with_named_internal_nodes(self):
        tc  = DendropyTreeCleaner()
        tree = Tree.get(data="('Asulf_Archaeoglobus.1_2280~2522125074':7.17,(('Afulgi_764~2528311132':0.0,'CP006577_764~2588253768':0.0):0.0,'AE000782_746~638154502':0.0)'s__Archaeoglobus fulgidus':7.555):1.461;\n",
                        schema='newick')
        tc.remove_sequences(tree,
                            ['CP006577_764~2588253768',
                             'Afulgi_764~2528311132'])
        self.assertEqual("(Asulf_Archaeoglobus.1_2280~2522125074:7.17,AE000782_746~638154502:7.555):1.461",
                         str(tree))

if __name__ == "__main__":
    unittest.main()
