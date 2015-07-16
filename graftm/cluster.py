
# Requires usearch/8.0.1623
import os
import logging
from collections import Counter
from graftm.sequence_io import SequenceIO

class readsClusterer:
    
    def __init__(self): 
        self.seqio = SequenceIO()
    
    def cluster(self, input_reads_path):
        '''
        cluster - input reads (unaligned) and output clustered fasta (unaligned)
        
        Parameters
        ----------
        input_reads_path : str
            path to input reads that are to be clustered (fasta)
        
        Returns
        -------
        output_reads_path : str
            path to clustered reads (fasta)
        '''
        
        nonredundant_sequences = {}
        seen=[]
        output=''
        
        logging.info('Clustering reads')
        reads=self.seqio.read_fasta_file(input_reads_path)
        logging.info('Read in %i reads' % len(reads))
        
        for read in reads:
            if read.seq in nonredundant_sequences:
                nonredundant_sequences[read.seq].append(read.name)
            else:
                nonredundant_sequences[read.seq] = [read.name]
                output+='%s\n%s\n' % (read.name, read.seq)
        
        
        
        path, basename    = os.path.split(input_reads_path)
        basename          = basename.split('.')[0]
        output_reads_path = os.path.join(path, basename+'_clustered.fa')
        
        with open(output_reads_path, 'w') as out:
            out.write(output)
        
        return output_reads_path, nonredundant_sequences
