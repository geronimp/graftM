import skbio.tree

class TreeCleaner:
    def clean_newick_file_for_fasttree_input(self, tree):
        '''Modify the given skbio.tree.TreeNode so that it can be output and fed
        into FastTree'''
        for n in tree.non_tips(include_self=True): n.name=None
        #convert underscores to spaces so skbio converts them back
        # without using single quotes for the names. See
        # https://github.com/biocore/scikit-bio/issues/934
        for n in tree.tips(): n.name = n.name.replace('_', ' ')
            
    def match_alignment_and_tree_sequence_ids(self, sequence_names, tree):
        '''Check to make sure that the sequences specified in the alignment
        and the tree are the same, otherwise raise an Exception detailing
        the problem for the user to fix'''
        
        tip_names_count = {}
        for t in tree.tips():
            # replace spaces with underscores as skbio interprets unquoted underscores
            # as spaces, where fasttree doesn't (I think)
            name = t.name.replace(' ','_')
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
        '''Remove sequences with in the given sequence_names array from the tree
        in place.
        
        Assumes the sequences are found in the tree, and that they are all unique.
        
        Potentially this could be better implemented using TreeNode.shear,
        but there are underscores to deal with (or not)
        '''
        for s in sequence_names:
            try:
                n = tree.find(s.replace('_',' '))
            except skbio.tree.MissingNodeError:
                # Just delete it already
                n = tree.find(s)
            if not n.parent.remove(n):
                raise Exception("Unexpectedly failed to remove '%s' from the\
                    tree. As a guess, are all the tip names unique?")
            tree.prune() #call prune() afer each remove as workaround for
            # https://github.com/biocore/scikit-bio/issues/970
        
        
        
    