from graftm.deduplicator import Deduplicator
from graftm.sequence_io import SequenceIO
from graftm.orfm import OrfM
import logging
import os

class Clusterer:

    def __init__(self):
        self.clust = Deduplicator()
        self.seqio = SequenceIO()
        self.seq_library = {}

        self.orfm_regex = OrfM.regular_expression()

    def uncluster_annotations(self, input_annotations, reverse_pipe):
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

        Returns
        -------
        output_annotations : hash
            An updated version of the above, which includes all reads from
            each cluster
        '''
        output_annotations = {}
        for placed_alignment_file_path, clusters in self.seq_library.iteritems():

            if reverse_pipe and placed_alignment_file_path.endswith("_reverse_clustered.fa"): continue
            placed_alignment_file = os.path.basename(placed_alignment_file_path)
            cluster_classifications = input_annotations[placed_alignment_file]

            if reverse_pipe:
                placed_alignment_base = placed_alignment_file.replace('_forward_clustered.fa', '')
            else:
                placed_alignment_base = placed_alignment_file.replace('_clustered.fa', '')
            output_annotations[placed_alignment_base] = {}
            for rep_read_name, rep_read_taxonomy in cluster_classifications.iteritems():

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
            output_path  = input_fasta.replace('_hits.aln.fa', '_clustered.fa')
            cluster_dict = {}

            logging.debug('Clustering reads')
            if os.path.exists(input_fasta):
                reads=self.seqio.read_fasta_file(input_fasta) # Read in FASTA records
                logging.debug('Found %i reads' % len(reads)) # Report number found
                clusters=self.clust.deduplicate(reads) # Cluster redundant sequences
                logging.debug('Clustered to %s groups' % len(clusters)) # Report number of clusters
                logging.debug('Writing representative sequences of each cluster to: %s' % output_path) # Report the name of the file
            else:
                logging.debug("Found no reads to be clustered")
                clusters = []

            self.seqio.write_fasta_file(
                                        [x[0] for x in clusters],
                                        output_path
                                        ) # Choose the first sequence to write to file as representative (all the same anyway)
            for cluster in clusters:
                cluster_dict[cluster[0].name]=cluster # assign the cluster to the dictionary
            self.seq_library[output_path]= cluster_dict

            output_fasta_list.append(output_path)

        return output_fasta_list


