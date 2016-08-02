from dendropy import Tree

class DendropyTreeCleaner:
    def write_fasttree_newick(self, tree, output_file_io):
        '''Write the given dendropy tree to the open output file I/O, ensuring that is
        acceptable for FastTree by e.g. removing quoted taxon names.

        Parameters
        ----------
        tree: dendropy.Tree
            tree to be written
        output_file_io: stream
            stream to write the newick formatted tree to

        '''
        newick = tree.as_string(schema='newick',
                                suppress_internal_node_labels=True,
                                suppress_internal_taxon_labels=True,
                                preserve_spaces=False) # Labels must not be quoted.
        output_file_io.write(newick.replace("'",'').replace(' ','_'))
