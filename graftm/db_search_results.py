from graftm.deduplicator import Deduplicator
from graftm.sequence_io import SequenceIO

import logging
import os

class DBSearchResult:
    '''Class for containing results from search pipeline in GraftM '''
    def __init__(self, ouput_reads, search_result, hit_read_count, slash_endings):
        self.clust = Deduplicator()
        self.seqio = SequenceIO()
        
        self.ouput_reads   = ouput_reads
        self.search_result = search_result
        self.hit_count     = hit_read_count
        self.slash_endings = slash_endings
    
    def hit_fasta(self):
        return self.ouput_reads
    
    def search_objects(self):
        return self.search_result
    
    def hit_count(self):
        return self.hit_count  
    
    def slash_endings(self):
        return self.slash_endings
    
    def _write_cluster(self, clusters):
        clustered_reads_output_path=self.ouput_reads.replace('.fa','_clustered.fa')
        logging.debug('Writing representative sequences of each cluster to: %s' % clustered_reads_output_path)
        cluster_dict={}
        with open(clustered_reads_output_path, 'w') as out:
            for cluster in clusters: 
                out.write( ">%s\n%s\n" % (cluster[0].name, cluster[0].seq) )
                cluster_dict[cluster[0].name]=cluster
        return clustered_reads_output_path, cluster_dict
    
    def cluster(self):
        logging.debug('Clustering reads')
        reads=self.seqio.read_fasta_file(self.hit_fasta())
        logging.debug('Found %i reads' % len(reads))
        self.groups=self.clust.deduplicate(reads)
        logging.debug('Clustered to %s groups' % len(self.groups))
        self.ouput_reads, self.clusters=self._write_cluster(self.groups)
        
    