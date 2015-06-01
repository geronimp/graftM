from skbio import TreeNode

class TreeCleaner:
    def clean_newick_file(self, input_tree_filename, output_tree_filename):
        tree = TreeNode.from_newick(open(input_tree_filename).read())
        with open(output_tree_filename,'w') as f:
            f.write(tree.to_newick())