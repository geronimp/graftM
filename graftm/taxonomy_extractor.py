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
                if n.name and n.parent:
                    tax = [n.name]+tax
                n = n.parent
            tip_name = tip.name.replace(' ','_')
            if tip_name in tip_to_taxonomy:
                raise Exception("Found duplicate tip name '%s'" % tip_name)
            if len(tax)==0:
                tip_to_taxonomy[tip_name] = []
            else:
                tip_to_taxonomy[tip_name] = [t.strip() for t in '; '.join(tax).split(';')]
        return tip_to_taxonomy