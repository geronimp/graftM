#!/usr/bin/env python

from skbio import TreeNode

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2014"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd, Ben Woodcroft"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"


class Rerooter:


    def importer(self, tree_path):
        t = TreeNode.from_newick(open(tree_path, 'r'))
        return t

    def node_list(self, t, nodes_path):
        id_list = []
        misses = []
        for i in open(nodes_path).readlines():
            try:
                id_list.append(t.find(i.rstrip()))
            except:
                misses.append(i.rstrip())

        if len(id_list) < 1:
            print 'No nodes to continue with. Exiting'
            exit()
        if len(misses) > 0:
            print 'The following nodes were not found in the tree, continuing with %s nodes:\n%s' % (len(id_list), str(' '.join(misses)))


        return id_list

    def lca(self, t, n):
        lca = t.lowest_common_ancestor(n)
        return lca

    def reroot(self, lca, t):
        root = lca.root()
        return root

    def write_tree(self, t, o):
        t.write(o, format='newick')

    def main(self, i, o, n):

        #load tree
        tree = self.importer(i)

        # Get nodes
        nodes = self.node_list(tree, n)

        # Find last common ancestor
        lca = self.lca(tree, nodes)

        # Reroot dat tree
        rerooted_tree = self.reroot(lca, tree)

        # Write to file
        self.write_tree(rerooted_tree , o)


