

class Deduplicator:
    r"""Deduplicates sequences"""
    
    def deduplicate(self, aligned_sequence_objects):
        r"""Sort the given aligned_sequence objects into an array of arrays,
        where input sequences are grouped iff they have the same sequence
        
        Parameters
        ----------
        aligned_sequence_objects: array of Sequence objects
            input sequences
            
        Returns
        -------
        Array of arrays of Sequence objects"""
        
        sequence_to_groups = {}
        for s in aligned_sequence_objects:
            try:
                sequence_to_groups[s.seq].append(s)
            except KeyError:
                sequence_to_groups[s.seq] = [s]
        return sequence_to_groups.values()
    
    def lca_taxonomy(self, deduplicated_sequences, taxonomy_hash):
        r'''Given a set of deduplicated sequences and a taxonomy hash,
        return the respective LCAs of taxonomy
        
        Parameters
        ----------
        deduplicated_sequences: Array of arrays of Sequence objects
            as output from deduplicate()
        taxonomy_hash: dictionary 
            of sequence names to taxonomy array (i.e. array of str)
        
        Returns
        -------
        Array of taxonomy LCAs'''
        
        to_return = []
        for dup_group in deduplicated_sequences:
            lca = taxonomy_hash[dup_group[0].name]
            for s in dup_group[1:]:
                for i, tax in enumerate(taxonomy_hash[s.name]):
                    if i >= len(lca) or tax != lca[i]:
                        lca = lca[:i]
                        break
                if len(lca) > len(taxonomy_hash[s.name]):
                    lca = lca[:len(taxonomy_hash[s.name])]
            to_return.append(lca)
        return to_return
                
            
        
        
        
        
