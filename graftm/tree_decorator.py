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
import tempfile
import subprocess
import numpy as np
import os
import random
from scipy import stats
from itertools import combinations
from collections import Counter
from skbio.tree import TreeNode
from graftm.taxonomy import Taxonomy

# Local imports
################################################################################
################################################################################
################################################################################

class TreeUnrootedException(Exception): 
    pass

#-------------------------------------------------------------------------------

class MalformedTreeException(Exception):
    pass

#-------------------------------------------------------------------------------

class TreeDecorator: # Room to grow!
    '''
    Decorate trees with information!
    '''
    GREENGENES_PREFIXES = set(['k__',
                               'p__',
                               'c__',
                               'o__',
                               'f__',
                               'g__',
                               's__'])

    def __init__(self, tree, no_unique_names):
        '''
        Parameters
        ----------
        tree: skbio.tree.TreeNode
            Tree to be decorated
        no_unique_names : bool
            True of False, whether names will be made non-redundant if 
            monophyletic groupings are split.
        '''
        self.tree               = tree
        self.no_unique_names    = no_unique_names
    
    def write_taxonomy(self, output_taxonomy_path):
        '''
        Write a greengenes formatted taxonomy file to the "output_taxonomy_path"
        path for each tip in the tree. If no classification is possible the 
        string will be left empty
        
        Parameters
        ----------
        output_taxonomy_path
        '''
        logging.debug("Writing decorated taxonomy to %s" % output_taxonomy_path)
        decorated_taxonomy = {}
        for tip in self.tree.tips():
            taxonomy_string = []
            tip_name = tip.name.replace(' ','_')
            index = 0
            for ancestor in tip.ancestors():
                if hasattr(ancestor, 'index'):
                    index += ancestor.index
                    taxonomy_string.append(ancestor.name)
                    if index > self.max_taxonomy_index:break
            taxonomy_string = reversed(taxonomy_string)
            decorated_taxonomy[tip_name] = Taxonomy.GG_SEP_0\
                                                   .join(taxonomy_string)
        
        with open(output_taxonomy_path, 'w') as output_taxonomy_io:
            for taxa_name, taxonomy_string in decorated_taxonomy.iteritems():
                line = "%s\t%s\n" % (taxa_name, taxonomy_string)
                output_taxonomy_io.write(line)
    
    def taxonomy_partition(self, taxonomy):
        '''
        Partition the tree into strict taxonomy clades. This code currently 
        does not allow for inconsistency within a clade of any fashion, but 
        should be flexible enough to integrate easily should the need arise.
        Taxonomy decoration should be granular, but if there are many HGT 
        events or the tree just doesn't follow taxonomy it could be messy. 
        A simple tally system distinguished taxonomic groupings. This can be 
        turned off by providing a false to the no_unique_names paramters in the
        __init__ of this function
        
        Parameters
        ----------

        
        taxonomy: Taxonomy()
                 
        Returns
        -------
        skbio TreeNode object
        '''
        self.max_taxonomy_index = taxonomy.max
        taxonomy_to_next_unique = {}
        to_decorate = []
        
        for node in self.tree.preorder():    
            tax_string_array = []
            
            for rank in range(self.max_taxonomy_index):
                node_classification = None
                consistent_annotation = True
                for tip in node.tips():
                    tip_name = tip.name.replace(' ', '_')
                    if tip_name in taxonomy.taxonomy:
                        classification = taxonomy.taxonomy[tip_name][rank]
                        if node_classification!=None:
                            if node_classification!=classification:
                                consistent_annotation=False
                        else:
                            node_classification=classification
                
                if consistent_annotation:
                    if self.no_unique_names:
                        tax_string_array.append(classification)
                    else:
                        if classification in taxonomy_to_next_unique:
                            unique_number = taxonomy_to_next_unique[classification]
                            tax_string_array.append('%s_%i' % (classification,
                                                               unique_number))
                            taxonomy_to_next_unique[classification]+=1
                        else:
                            tax_string_array.append(classification) 
                            taxonomy_to_next_unique[classification] = 1

            # At this point we have the have the list of classifications 
            # assigned to the node we are on e.g. ['k__Bacteria', 'p__Proteo...]
            
            # Here I'm stripping out any empty gg prefixes, if they exist.
            tax_string_array = [tax for tax in tax_string_array 
                                if tax not in self.GREENGENES_PREFIXES]
            
            if any(tax_string_array):
                # We need to keep track of the index we're on (i.e. which
                # rank we are up to in the taxonomy string.) so that we do not
                # Exceed the maximum number of ranks expected. For example if 
                # You had a tree where some bacteria are within an archaeal 
                # clade, the full tax string will not be joined to that of the
                # clase that it is within (i.e. the tax string will not be:
                # [k__Archaea, k__Bacteria, p__Proteo...]). Above we account 
                # For this situation by not allowing for any inconsistent 
                # clades but having this index allows us to implement a 
                # threshold, rather than strict clades, somewhere down the line
                index = 0  
                
                for anc in node.ancestors():
                    
                    # Here we count the number of indexes we encounter in the 
                    # ancestors of the node. 
                    if hasattr(anc, 'index'):
                        # If we exceed the number we expect
                        if index >= self.max_taxonomy_index:break
               
                # Then we strip off those extra taxonomy ranks (the k__Archaea 
                # in [k__Archaea, k__Bacteria, p__Proteo...])
                tax_string_array = tax_string_array[index:]
                # And if we have any tax string left
                if len(tax_string_array) > 0:
                    # Decorate the node with that information
                    node.name =  Taxonomy.GG_SEP_0.join(tax_string_array) 
                    node.index = len(tax_string_array)
        
    def return_tree(self):
        '''
        Returns the tree after all manipulation is done
        '''
        return self.tree                   
                    
################################################################################d
################################################################################
        