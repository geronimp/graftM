################################################################################
################################# - Imports - ##################################

# Local imports
from graftm.getaxnseq import Getaxnseq
from graftm.rerooter import Rerooter
from graftm.greengenes_taxonomy import GreenGenesTaxonomy
# System imports
from skbio import TreeNode
import argparse
import logging
import os
import sys

################################################################################
################################ - Statics - ###################################
MAX_RESOLUTION = 7

################################################################################
################################## - Code - ####################################

class TreeUnrootedException(Exception): 
    pass

class TreeDecorator:
    '''Yet another class that decorates trees, only the way I want it.'''
    
    
    def __init__(self, 
                 tree,
                 taxonomy,
                 seqinfo=None):
        '''
        Set up empty lists and open and import tree object, and load the
        taxonomy into a hash
        '''
        
        self.gg_prefixes = ["k__", "d__", 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        
        # Set empty list to record nodes that have already been assigned
        # taxonomy.
        self.encountered_nodes = {}
        self.encountered_taxonomies = set()
        
        # Read in tree
        self.tree = tree
        
        # Read in taxonomy
        
        if seqinfo:
            gtns = Getaxnseq()
            logging.info("Importing taxtastic taxonomy from files: %s and %s" % (taxonomy, seqinfo))
            self.taxonomy =  gtns.read_taxtastic_taxonomy_and_seqinfo(open(taxonomy), open(seqinfo))
            
            for id, taxonomy_list in self.taxonomy.iteritems():
                if len(taxonomy_list) != MAX_RESOLUTION:
                    self.taxonomy[id] = taxonomy_list + self.gg_prefixes[1:][len(taxonomy_list):]
        else:
            logging.info("Importing greengenes taxonomy from file: %s" % (taxonomy))
            self.taxonomy = GreenGenesTaxonomy.read_file(taxonomy).taxonomy
            
    def _nodes(self):
        '''
        Iterate through nodes in tree, returning those that are not tips (i.e.
        returning only nodes)
        
        Inputs:
            None
        Yields:
            Node object, if it is not a tip
        '''
        logging.debug("Iterating through nodes")       
        if len(self.tree.children) > 2:
            logging.warning("There are > 2 descendants from the root of the \
tree that was  provided. Tree is being rerooted to the branch of the node that \
is the greatest distance to the root.")
            self.tree = Rerooter().reroot(self.tree)
        for node in self.tree.preorder():
            if node.is_root():
                logging.debug("%i reads in tree from root" % (len(list(node.tips()))))                   
                yield node
            elif node.is_tip():
                logging.debug("Tip %s reached" % (node.name))
            else:
                if node not in self.encountered_nodes: # If node hath not been decorated.
                    if node.name:
                        logging.debug("Node %s reached" % (node.name))
                    else:
                        logging.debug("Unnamed node reached")
                    yield node
    
    def _transpose_taxonomy(self, taxonomy_list):
        '''
        Transpose the taxonomy grid to contain lists of each taxonomic
        rank
        
        Inputs: 
            taxonomy_list: array
                - List of split taxonomy strings to be tansposed
        Outputs: 
            transposed_taxonomy: array
                - As per above description
        '''
        transposed_taxonomy = []
        for index in range(0, 7):

            rank_taxonomy_list = [taxonomy_entry[index] 
                                    for taxonomy_entry in taxonomy_list]
            transposed_taxonomy.append(set(rank_taxonomy_list))
        return transposed_taxonomy
    
    def _write_tree(self, output):
        '''Just Writes tree to file uisng scikit-bio'''
        logging.info("Writing decorated tree to file: %s" % (output))
        self.tree.write(
                    output,
                    format = "newick"   
                        )
        
    def _write_consensus_strings(self, output):
        '''Writes the taxonomy of each leaf to a file'''
        logging.info("Writing decorated taxonomy to file: %s" % (output))
        with open(output, 'w') as out:
            for tip in self.tree.tips():
                ancestors = tip.ancestors()
                
                ancestor_tax = []
                for ancestor in ancestors:
                    if ancestor.name:
                        split_node_name = ancestor.name.split(':')
                        if len(split_node_name) == 2:
                            ancestor_tax.append(split_node_name[1])
                        elif len(split_node_name) == 1:
                            pass
                        else:
                            raise Exception("Malformed node name: %s" % ancestor.name)
                        
                tax_list = list(reversed(ancestor_tax))
                tax_name = tip.name.replace(" ", "_")
                if len(tax_list) < 1:
                    logging.debug("Species %s not decorated within the tree, using provided tax string as substitute" % (tax_name))
                    if tax_name in self.taxonomy:
                        tax_string = '; '.join(self.taxonomy[tax_name])
                    else:
                        logging.warning("No taxonomy found for species %s!" % (tax_name))
                        tax_string = "Unknown"
                else:
                    tax_string = '; '.join(tax_list)
                output_line = "%s\t%s\n" % (tax_name, tax_string)
                out.write(output_line)
    
    def _get_tax_index(self, ancestry):
        '''
        Iterate through ancestor nodes to find the current taxonomy index (i.e.
        the rank that should be assigned to the current node.
        
        Inputs:
            ancestry: list
                TreeNode.ancestory list that contains all nodes above the 
                current.
        Outputs:
            Current index of taxonomy. Returns None if node has no ancestors 
            (i.e. is root)
                
        '''
        current_index = 0
        ancetor_rank_number = 0
        ancestry = [x for x in ancestry if not x.is_root()]
        if any(ancestry):
            ancestor_tax_list = []
            for ancestor in ancestry:
                ai = self.encountered_nodes[ancestor]
                
                if ancestor.name:
                    split_node_name = ancestor.name.split(':')
                    if len(split_node_name) == 2:
                        ancestor_tax_list.append(split_node_name[1])
                    elif len(split_node_name) == 1:
                        pass
                    else:
                        raise Exception("Malformed node name: %s" % ancestor.name)
                if ai > current_index:
                    current_index = ai
            
            ancestor_tax = '; '.join(reversed(ancestor_tax_list))
            logging.debug("Ancestors at %i nodes: %s" % (len(ancestor_tax_list), ancestor_tax))
            ancetor_rank_number = len(ancestor_tax.split('; ')) #
            return current_index, ancetor_rank_number
        else:       
            return current_index, ancetor_rank_number
            #raise Exception("Programming error in _get_tax_index. Failed to find ancestors for current node.")
    
    def extract(self, output_tax):
        self._write_consensus_strings(output_tax)
    
    def decorate(self, output_tree, output_tax, no_unique_names):
        '''
        Main function for TreeDecorator class. Accepts the output directory, 
        iterates through nodes, classifying them according to provided taxonomy. 
        '''
        # Define list of prefixes        
        logging.info("Decorating tree")
        for node in self._nodes():
            placement_depth = 0 
            current_index = 0
            parent = node.parent
            children = list(node.tips())
            logging.debug("Current node has %i children" % (len(children)))
                        
            # Gather the children 
            children_taxonomy = []
            for t in children:
                tip_name = t.name.replace(' ', '_')
                if tip_name in self.taxonomy:
                    children_taxonomy.append(self.taxonomy[tip_name])
            node_tax = self._transpose_taxonomy(children_taxonomy)
            tax_list = []
            for index, rank in enumerate(node_tax):
                if len(rank) == 1:
                    tax=list(rank)[0]
                    if tax not in self.gg_prefixes:
                        tax_list.append(list(rank)[0])
                        placement_depth += 1
                else:
                    break

            if len(tax_list) >= 1:
                logging.debug("Node has consistent taxonomy: %s" % '; '.join(tax_list))
                
                if not parent:
                    current_index = len(tax_list)
                    resolution = 1
                else:
                    if parent.is_root():
                        if parent.name:
                            current_index, resolution = self._get_tax_index(node.ancestors())
                        else:
                            current_index = 0
                            resolution = len(tax_list)   
                    else:
                        current_index, resolution = self._get_tax_index(node.ancestors())
                        
                    if len(tax_list) <= current_index:
                        logging.debug("No consistent taxonomy found for node! Moving to next node.")
                        self.encountered_nodes[node] = current_index
                        continue
                    else:
                        if len(tax_list) == placement_depth:

                            tax_list = tax_list[current_index:]
                            for child in children:
                                self.encountered_nodes[child] = current_index
                        elif len(node.children) == len(list(node.tips())):
                            tax_list = tax_list[current_index:]
                            for child in children:
                                self.encountered_nodes[child] = current_index
                        else:
                            current_index = len(tax_list)
                
                if resolution < MAX_RESOLUTION:
                    new_tax_list = [] 
                    
                        
                    for tax in tax_list:
                        if not no_unique_names:
                            if tax in self.encountered_taxonomies:
                                idx = 1
                                while tax in self.encountered_taxonomies:
                                    if idx > 1: 
                                        tax = '%s_%i' % ('_'.join(tax.split('_')[:-1]), idx)
                                    else:
                                        tax = tax + '_%i' % idx
                                    idx += 1

                        new_tax_list.append(tax)
                        self.encountered_taxonomies.add(tax)
                    tax_string = '; '.join(new_tax_list)
                    
                    if node.name:
                        split_node_name = node.name.split(':')
                        if len(split_node_name) > 2: 
                            raise Exception("Malformed tree. Please remove ':' from the read names in the tree")
                            exit(1)
                        elif len(split_node_name) < 2:
                            node_name = "%s:%s" % (split_node_name[0], tax_string)
                        else:
                            bootstrap_value, node_name = split_node_name[0], split_node_name[1]
                            node_name = "%s:%s" % (bootstrap_value, tax_string)
                    else:
                        node_name = tax_string
                    node.name = node_name
                    logging.debug("Renamed node %s" % node_name)
                    self.encountered_nodes[node]=index
                    
                else:
                    self.encountered_nodes[node]=current_index
            else:
                logging.debug("Cannot resolve node, no consistent taxonomy beyond that which has been described.")
                self.encountered_nodes[node]=current_index
            
        if output_tree:
            self._write_tree(output_tree)
        if output_tax:
            self._write_consensus_strings(output_tax)

################################################################################
################################################################################
        
        