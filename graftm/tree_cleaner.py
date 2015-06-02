from skbio import TreeNode
import re

class TreeCleaner:
    def clean_newick_file(self, input_tree_filename, output_tree_filename):
        tree = TreeNode.read(open(input_tree_filename))
        for n in tree.non_tips(include_self=True): n.name=None
        r = re.compile(' ')
        with open(output_tree_filename,'w') as f:
            # the re.sub is required because skbio removes underscores
            # from names even if they were in the original newick file.
            f.write(r.sub('_', tree.to_newick(escape_name=False)))