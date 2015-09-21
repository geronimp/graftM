
# System imports
import logging
from skbio import TreeNode
from graftm.reannotator import Reannotator
from graftm.tree_decorator import TreeDecorator

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
            retroot the tree tree provided to old_tree
        old_tree_path: str
            Path to the file containing the tree to be re-rooted. This tree will
            be rerooted at the same position as the tree porovided to the 
            reference_tree
        '''
        reference_tree_path = kwargs.pop('reference_tree_path', None)
        tree_path = kwargs.pop('tree_path', None)
        
        if tree_path == None:
            logging.error("No tree for provided for Decorator class to modify")
            exit(1)
        else:
            logging.debug("Importing old tree from file: %s" 
                            % tree_path)
            self.tree = TreeNode.read(open(tree_path, "r"))  
        
        if reference_tree_path:
            logging.debug("Importing reference tree from file: %s" 
                            % reference_tree_path)
            self.reference_tree = TreeNode.read(open(reference_tree_path, "r"))   
        else:
            self.reference_tree = None
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments provided to Decorator class: %s" % kwargs)
        
    def _reroot(self):
        '''Run the re-rooting algorithm in the Reannotator class.'''
        reannotator = Reannotator()
        self.tree = reannotator._reroot_tree_by_old_root(self.reference_tree, 
                                                         self.tree)
        
    def main(self, taxonomy, output_tree, output_tax, seqinfo = None):
        '''Decorate and if necessary, re-root the tree. If an old reference tree
        is provided it is assumed that re-rooting is desired
        
        Parameters
        ----------
        taxonomy: str
            Path to a GreenGenes formatted taxonomy file, containing the 
            taxonomy of all or some if the sequences used to construct the tree.
            This taxonomy will be used by the tree_decorator.TreeDecorator class
            to decorate self.tree
        output_tree: str
            Path to file to which the decorated tree will be written to.
        output_tax: str
            Path to file to which the decorated taxonomy will be written to.'''        
        
        # Reroot
        if self.reference_tree:
            self._reroot()
            
        # Decorate
        td =  TreeDecorator(self.tree,
                            taxonomy,
                            seqinfo)
        td.decorate(output_tree, output_tax)
        