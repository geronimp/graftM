from graftm.rerooter import Rerooter

class TreeParaphyleticException(Exception): pass

class Reannotator:
    
    def _reroot_tree_by_old_root(self, old_tree, new_tree):
        r = Rerooter().reroot_by_tree(old_tree, new_tree)
        return r