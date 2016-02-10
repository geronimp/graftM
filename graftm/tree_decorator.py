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
__author__ = "Joel Boyd, Ben Woodcroft"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
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
        tree: skbio TreeNode obj
            http://scikit-bio.org/docs/latest/generated/skbio.tree.TreeNode.html#skbio.tree.TreeNode
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

        
        taxonomy: TaxonomyObject
                 
        Returns
        -------
        skbio TreeNode object
        '''
        self.max_taxonomy_index=taxonomy.max
        consistency_check = 1
        seen_taxonomy_rank = {}
        to_decorate = []
        
        for node in self.tree.preorder():    
            tax_string_array = []
            for rank in range(0, self.max_taxonomy_index):
                node_classifications = []
                for tip in node.tips():
                    tip_name = tip.name.replace(' ','_')
                    if taxonomy.check(tip_name):
                        node_classifications.append(taxonomy.rank(tip_name, 
                                                                  rank))
                    else:
                        to_decorate.append(node_classifications)

                if len(set(node_classifications)) == consistency_check:
                    classification = node_classifications[0]
                    if self.no_unique_names:
                        tax_string_array.append(classification)
                    else:
                        if classification in seen_taxonomy_rank:
                            unique_number = seen_taxonomy_rank[classification]
                            tax_string_array.append('%s_%i' % (classification,
                                                               unique_number))
                            seen_taxonomy_rank[classification]+=1
                        else:
                            tax_string_array.append(classification) 
                            seen_taxonomy_rank[classification] = 1
                            
            tax_string_array = [tax for tax in tax_string_array
                                if tax not in self.GREENGENES_PREFIXES]
            
            if any(tax_string_array):
                index = 0
                for anc in node.ancestors():
                    
                    if hasattr(anc, 'index'):
                        index+=anc.index
                        if index >= self.max_taxonomy_index:break

                tax_string_array = tax_string_array[index:]
                
                if len(tax_string_array) >0:
                    node.name =  Taxonomy.GG_SEP_0.join(tax_string_array)
                    node.index = len(tax_string_array)
        
    def return_tree(self):
        '''
        Returns the tree after all manipulation is done
        '''
        return self.tree                   
                    
################################################################################d
################################################################################
        