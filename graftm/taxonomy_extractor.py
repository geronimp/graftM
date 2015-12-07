import re

class TaxonomyExtractor:
    def taxonomy_from_annotated_tree(self, tree_node):
        '''Given an annotated tree, return the taxonomy of each tip
        
        Parameters
        ---------
        tree_node: skbio.tree.TreeNode
        
        Returns
        -------
        dictionary of tip name to array of taxonomic affiliations'''
        tip_to_taxonomy = {}
        for tip in tree_node.tips():
            tax = []
            n = tip.parent
            while n:
                node_taxonomy = self.taxonomy_from_node_name(n.name)
                if node_taxonomy and n.parent:
                    tax = [node_taxonomy]+tax
                n = n.parent
            tip_name = tip.name.replace(' ','_')
            if tip_name in tip_to_taxonomy:
                raise Exception("Found duplicate tip name '%s'" % tip_name)
            if len(tax)==0:
                tip_to_taxonomy[tip_name] = []
            else:
                tip_to_taxonomy[tip_name] = [t.strip() for t in '; '.join(tax).split(';')]
        return tip_to_taxonomy
    
    def taxonomy_from_node_name(self, node_name):
        '''return the taxonomy incorporated at a particular node, or None
        if it does not encode any taxonomy
        
        Parameters
        ----------
        node_name: str
            a TreeNode.name
            
        Returns
        -------
        Taxonomy as a string, or None if there is no taxonomy
        '''
        
        def is_float(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
        
        if node_name is None:
            return None
        elif is_float(node_name):
            # no name, just a bootstrap
            return None
        else:
            bootstrap_regex = re.compile(r'[\d\.]+:(.*)')
            reg = bootstrap_regex.match(node_name)
            if reg:
                # bootstrap in name
                return reg.groups(0)[0]
            else:
                # bootstrap not in name
                return node_name
