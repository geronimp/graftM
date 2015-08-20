from graftm.deduplicator import Deduplicator
from graftm.sequence_io import SequenceIO

import logging

class DBSearchResult:
    '''
    DBSearchResult - Class for containing results from search pipeline in GraftM
    '''
    def __init__(self, output_reads, search_result, hit_read_count, slash_endings):
        self.clust = Deduplicator()
        self.seqio = SequenceIO()
        
        self.output_reads   = output_reads
        self.search_result = search_result
        self.hit_count     = hit_read_count
        self.slash_endings = slash_endings
    
    def hit_fasta(self): # Return the path to the fasta file of hits
        return self.output_reads
    
    def search_objects(self): # Return SequenceSearchResult object with parameters defined
        return self.search_result
    
    def hit_count(self): # Return the number of hits found
        return self.hit_count  
    
    def slash_endings(self): # Return True if sequence headers end in /1 or /2, False if not
        return self.slash_endings
    
    def _write_cluster(self, clusters): 
        '''
        _write_cluster - Writes representative sequences of each cluster to file
        and returns the name of the file written to.
        
        Parameters
        ----------
        clusters : list
            list of lists, one for each cluster containing the sequences within
            each cluster.
        
        Returns
        -------
        clustered_reads_output_path: str
            Path to the file containing representative sequences of each cluster
            in FASTA format.
            
        cluster_dict : hash
            A dictionary with the header of the representative sequence 
            pertaining to each cluster as the key, and a list of the sequence 
            objects (from SequenceIO with parameters defined) that belong to
            each cluster.
        '''
        clustered_reads_output_path=self.output_reads.replace('.fa','_clustered.fa') # Define the output path to write representative sequences to
        logging.debug('Writing representative sequences of each cluster to: %s' % clustered_reads_output_path) # Report the name of the file
        cluster_dict={} # Define an empty deictionary to store cluster information
        with open(clustered_reads_output_path, 'w') as out: 
            for cluster in clusters: # for every cluster
                out.write( ">%s\n%s\n" % (cluster[0].name, cluster[0].seq) ) # Choose the first sequence to write to file as representative (all the same anyway)
                cluster_dict[cluster[0].name]=cluster # assign the cluster to the dictionary
        return clustered_reads_output_path, cluster_dict 
    
    def cluster(self):
        '''
        cluster - Clusters reads at 100% identity level and  writes them to 
        file. Resets the hit_fasta variable as the FASTA file containing the 
        clusters.
        '''
        logging.debug('Clustering reads')
        reads=self.seqio.read_fasta_file(self.hit_fasta()) # Read in FASTA records
        logging.debug('Found %i reads' % len(reads)) # Report number found
        self.groups=self.clust.deduplicate(reads) # Cluster redundant sequences
        logging.debug('Clustered to %s groups' % len(self.groups)) # Report number of clusters
        self.output_reads, self.clusters=self._write_cluster(self.groups) # Write to file, and re-define output_reads and set cluster variable.
        
    