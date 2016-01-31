from graftm.deduplicator import Deduplicator
from graftm.sequence_io import SequenceIO
from graftm.orfm import OrfM
import logging
import os
import re

class Clusterer:

    def __init__(self, slash_endings):
        self.clust = Deduplicator()
        self.seqio = SequenceIO()
        self.seq_library = {}
        self.slash_endings = slash_endings
        self.orfm_regex = OrfM.regular_expression()


    def uncluster_annotations(self, input_annotations, reverse_pipe,
                              merge_reads):
        '''
        Update the annotations hash provided by pplacer to include all
        representatives within each cluster

        Parameters
        ----------
        input_annotations : hash
            Classifications for each representative sequence of the clusters.
            each key being the sequence name, and the entry being the taxonomy
            string as a list.
        reverse_pipe : bool
            True/False, whether the reverse reads pipeline is being followed.
        merge_reads : bool
            True/False, whether forward and reverse reads were merged for 
            placement

        Returns
        -------
        output_annotations : hash
            An updated version of the above, which includes all reads from
            each cluster
        '''
        output_annotations = {}
        for _, (clusters, placed_alignment_file_path) in self.seq_library.iteritems():

            if reverse_pipe and placed_alignment_file_path.endswith("_reverse_clustered.aln.fa"): continue
            placed_alignment_file = os.path.basename(placed_alignment_file_path)
            cluster_classifications = input_annotations[placed_alignment_file]

            if reverse_pipe:
                placed_alignment_base = placed_alignment_file.replace('_forward_clustered.aln.fa', '')
            else:
                placed_alignment_base = placed_alignment_file.replace('_clustered.aln.fa', '')
            output_annotations[placed_alignment_base] = {}
            for rep_read_name, rep_read_taxonomy in cluster_classifications.iteritems():
                if merge_reads:
                    rep_read_name=rep_read_name[:-2]
                if reverse_pipe:
                    orfm_regex = OrfM.regular_expression()
                    clusters={(orfm_regex.match(key).groups(0)[0] if orfm_regex.match(key) else key): item for key, item in clusters.iteritems()}
                for read in clusters[rep_read_name]:
                    output_annotations[placed_alignment_base][read.name] = rep_read_taxonomy

        return output_annotations

    def cluster(self, input_fasta_list, reverse_pipe):
        '''
        cluster - Clusters reads at 100% identity level and  writes them to
        file. Resets the input_fasta variable as the FASTA file containing the
        clusters.

        Parameters
        ----------
        input_fasta_list : list
            list of strings, each a path to input fasta files to be clustered.
        reverse_pipe : bool
            True/False, whether the reverse reads pipeline is being followed.

        Returns
        -------
        output_fasta_list : list
            list of strings, each a path to the output fasta file to which
            clusters were written to.
        '''
        
        output_fasta_list = []
        for input_fasta in input_fasta_list:
            output_path  = input_fasta.replace('_hits.aln.fa', '_clustered.aln.fa')
            path = os.path.split(input_fasta)[0]
            cluster_dict = {}

            logging.debug('Clustering reads')
            reads=self.seqio.read_fasta_file(input_fasta) # Read in FASTA records
            logging.debug('Found %i reads' % len(reads)) # Report number found
            clusters=self.clust.deduplicate(reads) # Cluster redundant sequences
            logging.debug('Clustered to %s groups' % len(clusters)) # Report number of clusters
            logging.debug('Writing representative sequences of each cluster to: %s' % output_path) # Report the name of the file

            self.seqio.write_fasta_file(
                                        [x[0] for x in clusters],
                                        output_path
                                        ) # Choose the first sequence to write to file as representative (all the same anyway)
            for cluster in clusters:
                if self.slash_endings:
                    name = cluster[0].name[:-2]
                else:
                    name = cluster[0].name                    
                cluster_dict[name]=cluster # assign the cluster to the dictionary
            self.seq_library[path]= [cluster_dict, output_path] 


            output_fasta_list.append(output_path)

        return output_fasta_list


