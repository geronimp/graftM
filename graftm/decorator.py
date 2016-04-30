#!/usr/bin/env python
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
# System imports
import logging
from skbio import TreeNode
from graftm.reannotator import Reannotator
from graftm.tree_decorator import TreeDecorator

# Local imports
################################################################################
################################################################################
################################################################################
class Decorator:
    '''re-root a tree and decorate it using old taxonomy, for the graftM 
    decorate pipeline.
    
    ALL trees provided to this class must be in newick format. No other formats 
    are supported. If any taxonomic decoration is already in the tree to be 
    decorated by this class is WILL be overwritten. Of course, bootstrap values 
    will remain untouched.'''

    def __init__(self, **kwargs):     
        '''
        Parameters
        ----------
        reference_tree_path: str
            Path to the file containing the reference tree, which is used to
            retroot the tree tree provided to tree
        old_tree_path: str
            Path to the file containing the tree to be re-rooted. This tree will
            be rerooted at the same position as the tree provided to the 
            reference_tree
        '''
        reference_tree_path = kwargs.pop('reference_tree_path', None)
        tree_path = kwargs.pop('tree_path')
        
        logging.debug("Importing old tree from file: %s" 
                        % tree_path)
        self.tree = TreeNode.read(open(tree_path))  
        
        if reference_tree_path:
            logging.debug("Importing reference tree from file: %s" 
                            % reference_tree_path)
            self.reference_tree = TreeNode.read(open(reference_tree_path, "r"))   
        else:
            self.reference_tree = reference_tree_path
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments provided to Decorator class: %s" % kwargs)
        
    def _reroot(self):
        '''Run the re-rooting algorithm in the Reannotator class.'''
        reannotator = Reannotator()
        self.tree = reannotator._reroot_tree_by_old_root(self.reference_tree, 
                                                         self.tree)
        
    def main(self, taxonomy, output_tree, output_tax, no_unique_tax, decorate):
        '''Decorate and if necessary, re-root the tree. If an old reference tree
        is provided it is assumed that re-rooting is desired
        
        Parameters
        ----------
        taxonomy: str
            Path to a GreenGenes formatted, or taxtastic formatted taxonomy 
            file, containing the taxonomy of all or some if the sequences used 
            to construct the tree. This taxonomy will be used by the 
            tree_decorator.TreeDecorator class to decorate self.tree
        output_tree: str
            Path to file to which the decorated tree will be written to.
        output_tax: str
            Path to file to which the decorated taxonomy will be written to.
            A boolean False if no output taxonomy file was specified. In this
            case no taxonomy is written to file.
        no_unique_tax: bool
            True or False indicating whether or not to make the annotations 
            assigned to each node within the tree non-redundant.
        decorate: bool
            Indicates that the tree is to be decorated with taxonomy.
        '''        
        
        # Reroot
        if self.reference_tree:
            self._reroot()
            
        # Decorate
        if decorate:
            td =  TreeDecorator(self.tree,
                                no_unique_tax)
            td.taxonomy_partition(taxonomy)
            
        if output_tax:
            td.write_taxonomy(output_tax)
        else:
            logging.info("Writing tree to file: %s" % output_tree)
            self.tree = td.return_tree()
            self.tree.write(
                            output_tree,
                            format = "newick"   
                            )