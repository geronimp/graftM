###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
#
#
#
###############################################################################
################################ - Imports - ##################################
# Local

# System
import itertools
import logging
import heapq
from dendropy import Tree

################################ - Statics - ##################################

################################ - Classes - ##################################

class BinaryTreeError(Exception):
    pass
class MalformedTreeError(Exception):
    pass
class TreeParaphyleticException(Exception):
    pass

################################# - Code - ####################################

class Rerooter:
    def __init__(self): pass
    
    def reroot(self, tree):
        '''
        Reroot a tree: This small piece of code is intended to take a rooted 
        tree that is NOT binary, that is to say, the root node has 3 or more
        child nodes. To create a binary tree, this function will find the two 
        longest distances between the child nodes, find which child node is 
        common to them both, and reroot at the branch connecting that node 
        to root. For example, the new root will be placed where the "x" is in 
        the following example tree:
                                      _
                                   __|_     _
                                  |_____x__|
                                --|__ _    |_
                                     |_
        
        Parameters
        ----------
        tree: dendropy.Tree
            Tree to be rerooted
            
        Raises
        ------
        BinaryTreeError : If the input tree is already a binary tree

        '''

        children = {x.length:x for x in tree.seed_node.child_edges()}
        tree.is_rooted=True
        
        if len(children) > 2:
            longest_edge_length = max(children.keys())
            longest_edge = children[longest_edge_length]
            tree.reroot_at_edge(longest_edge, 
                                length1=(longest_edge_length/2), 
                                length2=(longest_edge_length/2))
            return tree
        
        elif len(children) == 2:
            raise BinaryTreeError("Tree is already binary. Something has \
-mislead me here!")
        else:
            raise MalformedTreeError("Input tree is malformed")

    def reroot_by_tree(self, old_tree, new_tree):
        '''reroot the new tree so that it matches the old tree's root, if 
        possible. If more than one rerooting is possible, root at the longest
        internal branch that is consistent with the root of the old_tree.

        Assumes that the tree is binary. Both the old and new trees may be
        modified by this method.
        
        Parameters
        ----------
        old_tree: dendropy.Tree
            The old tree to try to match the root from
        new_tree: dendropy.Tree
            tree to be rerooted in the same place as the old_tree.
            Must include at least one leaf from each side of the old_tree's
            root (matching by node.name), but may not have all leaves from
            the old tree, and can have extra leaves. 
            
        Returns
        -------
        The new_tree rooted by the root of the old tree
        
        Exceptions
        ----------
        TreeParaphyleticException
            If either of the old_tree's branches are not monophyletic in the
            new tree

        '''
        
        # Ensure that the tree is rooted to avoid gotchas
        # e.g. https://github.com/jeetsukumaran/DendroPy/issues/51
        old_tree.is_rooted = True
        new_tree.is_rooted = True
        if len(old_tree.seed_node.child_nodes()) != 2:
            raise Exception("Unexpectedly found a non-binary tree. Perhaps need to use Rerooter.reroot() ?")
        # make a list of the left and right leaf names that are in the new tree
        new_tip_names = set([tip.taxon.label for tip in new_tree.leaf_node_iter()])
        old_left_tip_names = [tip.taxon.label for tip in old_tree.seed_node.child_nodes()[0].leaf_nodes() if tip.taxon.label in new_tip_names]
        old_right_tip_names = [tip.taxon.label for tip in old_tree.seed_node.child_nodes()[1].leaf_nodes() if tip.taxon.label in new_tip_names]
        
        # find the LCA of the lefts and the rights, giving up at the root.
        left_lca = new_tree.mrca(taxon_labels=old_left_tip_names)
        right_lca = new_tree.mrca(taxon_labels=old_right_tip_names)
        
        # if both didn't LCA before hitting the root, tree paraphyletic
        # take the first one where the LCA was hit, reroot here.
        # find the LCA of the other in the other half of the tree.
        # if failed before getting to the root, then tree paraphyletic
        # reroot on one side of the internal branch between the two LCAs
        if left_lca == new_tree.seed_node:
            if right_lca == new_tree.seed_node:
                raise TreeParaphyleticException("Tree paraphyletic case #1")
            else:
                new_tree.reroot_at_edge(right_lca.edge)
                new_lca = new_tree.mrca(taxon_labels=old_left_tip_names)
        else:
            new_tree.reroot_at_edge(left_lca.edge, length1=left_lca.edge.length)
            new_lca = new_tree.mrca(taxon_labels=old_right_tip_names)
            
        if new_lca.edge.rootedge:
            raise TreeParaphyleticException("Tree paraphyletic case #2")
        rerooting_edge = self._find_longest_internal_edge(new_lca)
        if rerooting_edge and rerooting_edge.head_node and rerooting_edge.tail_node:
            new_tree.reroot_at_edge(rerooting_edge, length1=rerooting_edge.length)
        return new_tree
        
    def _find_longest_internal_edge(self, node):
        '''return the node that has the longest branch length between the given
        node and the root
        
        Parameters
        ----------
        node: dendropy.Tree
            a node from the tree
            
        Returns
        -------
        The node that has the largest length between the node and the root_node
        '''
        max_length = -1
        max_edge = None
        while node and not node.edge.rootedge:
            if node.edge.length > max_length:
                max_edge = node.edge
                max_length = max_edge.length
            node = node.parent_node
        return max_edge
         
###############################################################################
###############################################################################
