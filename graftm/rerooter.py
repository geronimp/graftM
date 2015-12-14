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
from skbio import TreeNode


################################ - Statics - ##################################

################################ - Classes - ##################################

class BinaryTreeError(Exception):
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
        tree: TreeNode object
            TreeNode object (tree opened with skbio)
            
        Raises
        ------
        BinaryTreeError : If the input tree is already a binary tree

        '''   
        
        if len(tree.children)==2:
            raise BinaryTreeError("Tree is already binary. Something has \
-mislead me here!")
        
        length_dict = {node.length:node for node in tree.traverse()}
        tree = tree.root_at(length_dict[max(length_dict.keys())].parent)
        children = [child for child in tree.children if not child.is_tip()]
        no_found_root = True
        
        if len(children) > 2:
            
            node_combinations = list(itertools\
                                   .combinations([x for x in children], 
                                                 2))

            length_combination = list(itertools\
                                   .combinations([x.length for x in children], 
                                                 2))
            
            sums = [x[0] + x[1] for x in length_combination]
            
            while no_found_root:
                
                nodes = [node_combinations[sums.index(x)]
                         for x in heapq.nlargest(2, sums)]
                most_distant_node = [node for node in nodes[0] 
                                     if node in nodes[1]]
                if len(most_distant_node) > 0:
                    most_distant_node = most_distant_node[0]
                    no_found_root = False
                
            other_nodes = [node for node in children 
                           if node != most_distant_node]
            dec_len = most_distant_node.length/2
            most_distant_node.length = dec_len
            link_root = TreeNode(children=other_nodes, length = dec_len)
            tree = TreeNode(children=[link_root, most_distant_node], 
                            name = "root",
                            length=0)
            return tree
        
        else:
            logging.warning("This tree contains a root positioned \
at the base of one or more leaves. The root needs to be deeper for this type \
of re-rooting")
            logging.warning("This tree will be rooted at the mid-point as an \
an alternative.")
            tree=tree.root_at_midpoint()
            return tree
        
###############################################################################
###############################################################################
