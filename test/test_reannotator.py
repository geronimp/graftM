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

import os
import sys
import unittest
import logging
from skbio.tree._tree import TreeNode
from StringIO import StringIO

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from graftm.reannotator import Reannotator, TreeParaphyleticException

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    def assert_tree_equal_no_labels(self, expected_newick, observed_tree):
        '''should include some tree ordering because ordering of children is not relevant, but eh for now'''
        expected = TreeNode.read(StringIO(expected_newick))
        for n in expected.non_tips(): n.name = None
        for n in observed_tree.non_tips(): n.name = None
        self.assertEqual(str(expected), str(observed_tree))
        
    def test_hello_world(self):
        self.assert_tree_equal_no_labels('((C,(D,E):2.0):0.0,(A,B):4.0);\n',
             Reannotator()._reroot_tree_by_old_root(\
                TreeNode.read(StringIO('((A,B):1,(C,D):2);')),
                TreeNode.read(StringIO('((A,B):1,(C,(D,E):2):3);'))))
        
    def test_simple(self):
        self.assert_tree_equal_no_labels('((C:12.0,(D:13.0,E:14.0):2.0):0,(A:10.0,B:11.0):4.0);\n',
             Reannotator()._reroot_tree_by_old_root(\
                TreeNode.read(StringIO('((A,B):1,(C,D):2);')),
                TreeNode.read(StringIO('((C:12,(A:10,B:11)a:4)b:0.5,(D:13,E:14)c:1.5);'))))
        
    def test_extra_unannotated_at_root(self):
        self.assert_tree_equal_no_labels('((C:12.0,(D:13.0,E:14.0):2.0):4.0,(F:15.0,(A:10.0,B:11.0):1.0):0.0)root;\n',
             Reannotator()._reroot_tree_by_old_root(\
                TreeNode.read(StringIO('((A,B):1,(C,D):2);')),
                TreeNode.read(StringIO('((C:12,((A:10,B:11)d:1,F:15)a:4)b:0.5,(D:13,E:14)c:1.5);'))))
        
        self.assert_tree_equal_no_labels('((F:15.0,(C:12.0,a:4.0,(D:13.0,E:14.0):2.0)):0.0,(A:10.0,B:11.0):40.0);\n',
             Reannotator()._reroot_tree_by_old_root(\
                TreeNode.read(StringIO('((A,B):1,(C,D):2);')),
                TreeNode.read(StringIO('((C:12,((A:10,B:11)d:40,F:15),a:4)b:0.5,(D:13,E:14)c:1.5);'))))
        
    def test_paraphyletic_at_root(self):
        with self.assertRaises(TreeParaphyleticException):
            Reannotator()._reroot_tree_by_old_root(\
                TreeNode.read(StringIO('((A,B):1,(C,D):2);')),
                TreeNode.read(StringIO('((A,C):1,(C,B):2);')))



if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR)
    unittest.main()
    
