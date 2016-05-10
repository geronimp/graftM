import logging

from io import StringIO
from skbio.tree import TreeNode
from dendropy import Tree

from graftm.rerooter import Rerooter

class TreeParaphyleticException(Exception): pass

class Reannotator:
    
    def _reroot_tree_by_old_root(self, old_tree, new_tree):
        dendro_old = self.skbio_tree_to_dendropy(old_tree)
        dendro_new = self.skbio_tree_to_dendropy(new_tree)
        r = Rerooter().reroot_by_tree(dendro_old, dendro_new)
        return self.dendropy_tree_to_skbio(r)

    def skbio_tree_to_dendropy(self, skbio_tree):
        return Tree.get(
            data=str(skbio_tree),
            schema='newick')

    def dendropy_tree_to_skbio(self, dendropy_tree):
        return TreeNode.read(StringIO(unicode(str(dendropy_tree)+';')))
