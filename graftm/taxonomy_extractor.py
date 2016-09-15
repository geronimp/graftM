import re

class TaxonomyExtractor:
    def taxonomy_from_annotated_tree(self, dendropy_tree):
        '''Given an annotated dendropy tree, return the taxonomy of each tip
        
        Parameters
        ---------
        tree: dendropy.Tree
            Decorated tree to extract taxonomy from       
        Returns
        -------
        dictionary of tip name to array of taxonomic affiliations'''
        tip_to_taxonomy = {}
        for tip in dendropy_tree.leaf_node_iter():
            tax = []
            n = tip.parent_node
            while n:
                node_taxonomy = self.taxonomy_from_node_name(n.label)
                if node_taxonomy and n.parent_node:
                    tax = [node_taxonomy]+tax
                n = n.parent_node
            tip_name = tip.taxon.label.replace(' ','_')
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
            a node label.
            
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
