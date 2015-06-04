from skbio.tree import TreeNode
from Bio import SeqIO

class TreeCleaner:
    def clean_newick_file_for_fasttree_input(self, input_tree_filename, output_tree_filename):
        tree = TreeNode.read(open(input_tree_filename))
        for n in tree.non_tips(include_self=True): n.name=None
        #convert underscores to spaces so skbio converts them back
        # without using single quotes for the names. See
        # https://github.com/biocore/scikit-bio/issues/934
        for n in tree.tips(): n.name = n.name.replace('_', ' ')
        
        with open(output_tree_filename,'w') as f:
            tree.write(f)
            
    def match_alignment_and_tree_sequence_ids(self, alignment_file, newick_file):
        '''Check to make sure that the sequences specified in the alignment
        and the tree are the same, otherwise raise an Exception detailing
        the problem for the user to fix'''
        
        tip_names_count = {}
        tree = TreeNode.read(open(newick_file),format='newick')
        for t in tree.tips():
            # replace spaces with underscores as skbio interprets unquoted underscores
            # as spaces, where fasttree doesn't (I think)
            name = t.name.replace(' ','_')
            if name in tip_names_count:
                raise Exception("Duplicate tip name found in tree: '%s'" % name)
            else:
                tip_names_count[name] = 1
        for s in SeqIO.parse(open(alignment_file, 'r'), 'fasta'):
            name = s.name
            if name not in tip_names_count:
                raise Exception("The alignment sequence '%s' was found in the alignment but not the tree" % name)
            elif tip_names_count[name] > 1:
                raise Exception("Found duplicate sequence name '%s'" % name)
            else:
                tip_names_count[name] += 1
        for name, count in tip_names_count.iteritems():
            if count < 2:
                raise Exception("Sequence '%s' was found in the tree but not the alignment" % name)