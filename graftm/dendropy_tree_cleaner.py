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
            
    def match_alignment_and_tree_sequence_ids(self, sequence_names, tree):
        '''Check to make sure that the sequences specified in the alignment
        and the tree are the same, otherwise raise an Exception detailing
        the problem for the user to fix

        Parameters
        ----------
        sequence_names: list of str
            names of sequences to ensure are in the tree
        tree: dendropy.Tree
            tree to find names in

        '''
        
        tip_names_count = {}
        for t in tree.leaf_node_iter():
            # replace spaces with underscores as this is how they are given to FastTree.
            name = t.taxon.label.replace(' ','_')
            if name in tip_names_count:
                raise Exception("Duplicate tip name found in tree: '%s'" % name)
            else:
                tip_names_count[name] = 1
        for name in sequence_names:
            if name not in tip_names_count:
                raise Exception("The alignment sequence '%s' was found in the alignment but not the tree" % name)
            elif tip_names_count[name] > 1:
                raise Exception("Found duplicate sequence name '%s'" % name)
            else:
                tip_names_count[name] += 1
        for name, count in tip_names_count.iteritems():
            if count < 2:
                raise Exception("Sequence '%s' was found in the tree but not the alignment" % name)

    def remove_sequences(self, tree, sequence_names):
        '''Remove sequences with in the given sequence_names array from the tree in
        place. Assumes the sequences are found in the tree, and that they are
        all unique.

        Parameters
        ----------
        tree: dendropy.Tree
            tree to remove from
        sequence_names: list of str
            list of tip names to remove

        '''
        tree.prune_taxa_with_labels(sequence_names)
        tree.prune_taxa_with_labels([s.replace('_',' ') for s in sequence_names])
