from graftm.deduplicator import Deduplicator
from graftm.sequence_io import SequenceIO
import logging
import os

class Clusterer:
    
    def __init__(self):
        self.clust = Deduplicator()
        self.seqio = SequenceIO()
        self.seq_library = {}
    
    def uncluster_annotations(self, input_annotations):
        '''
        Update the annotations hash provided by pplacer to inlcude all 
        representatives within each cluster
        
        Parameters
        ----------
        input_annotations : hash
            Classifications for each representative sequence of the clusters. 
            each key being the sequence name, and the entry being the taxonomy
            string as a list. 
        
        Returns
        -------
        output_annotations : hash
            An updated version of the above, which includes all reads from 
            each cluster
        '''
        output_annotations = {}
        for placed_alignment_file_path, clusters in self.seq_library.iteritems():
            placed_alignment_file = os.path.basename(placed_alignment_file_path)
            cluster_classifications = input_annotations[placed_alignment_file]
            
            placed_alignment_base = placed_alignment_file.split('.')[0].replace('_hits', '')
            output_annotations[placed_alignment_base] = {}
            
            for rep_read_name, rep_read_taxonomy in cluster_classifications.iteritems():
                for read in clusters[rep_read_name]:
                    output_annotations[placed_alignment_base][read.name] = rep_read_taxonomy
        
        return output_annotations        
    
    def _write_cluster(self, clusters, output_path): 
        '''
        _write_cluster - Writes representative sequences of each cluster to file
        and returns the name of the file written to.
        
        Parameters
        ----------
        clusters : list
            list of lists, one for each cluster containing the sequences within
            each cluster.
        
        output_path: str
            Path to the file to which representative sequences of each cluster
            in FASTA format will be written.
        Returns
        -------
        cluster_dict : hash
            A dictionary with the header of the representative sequence 
            pertaining to each cluster as the key, and a list of the sequence 
            objects (from SequenceIO with parameters defined) that belong to
            each cluster.
        '''
        logging.debug('Writing representative sequences of each cluster to: %s' % output_path) # Report the name of the file
        cluster_dict={} # Define an empty deictionary to store cluster information
        with open(output_path, 'w') as out: 
            for cluster in clusters: # for every cluster
                out.write( ">%s\n%s\n" % (cluster[0].name, cluster[0].seq) ) # Choose the first sequence to write to file as representative (all the same anyway)
                cluster_dict[cluster[0].name]=cluster # assign the cluster to the dictionary
        self.seq_library[output_path]= cluster_dict
    
    def cluster(self, input_fasta_list):
        '''
        cluster - Clusters reads at 100% identity level and  writes them to 
        file. Resets the input_fasta variable as the FASTA file containing the 
        clusters.
        
        Parameters
        ----------
        input_fasta_list : list
            list of strings, each a path to input fasta files to be clustered.
        Retruns
        -------
        output_fasta_list : list
            list of strings, each a path to the output fasta file to which 
            clusters were written to.
        '''
        output_fasta_list = []
        for input_fasta in input_fasta_list:
            output_path = input_fasta.replace('.fa', '_clustered.fa')
            
            logging.debug('Clustering reads')
            reads=self.seqio.read_fasta_file(input_fasta) # Read in FASTA records
            logging.debug('Found %i reads' % len(reads)) # Report number found
            groups=self.clust.deduplicate(reads) # Cluster redundant sequences
            logging.debug('Clustered to %s groups' % len(groups)) # Report number of clusters
            clusters=self._write_cluster(groups, output_path) # Write to file, and re-define output_reads and set cluster variable.
            
            output_fasta_list.append(output_path)
        
        return output_fasta_list
        
        
        