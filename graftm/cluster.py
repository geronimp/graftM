
# Requires usearch/8.0.1623
import os
import subprocess
import logging

class readsClusterer:
    
    def __init__(self): pass   
    
    def _read_uc(self, uc):
        logging.debug('Reading in uc file: %s' % uc)
        uc_hash={}
        uc_lines=[line.strip().split() for line in open(uc).readlines()]
        for line in uc_lines:
            if line[0] == "S":
                ind=line[8]
                uc_hash[ind]=[]
            elif line[0] == "H":
                uc_hash[ind].append(line[8])
            elif line[0] == "C":
                continue
        logging.debug('%i unique clusters found' % len(uc_hash.keys()))
        return uc_hash
        
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
        logging.info('Clustering reads')
        path, basename    = os.path.split(input_reads_path)
        basename          = basename.split('.')[0]
        output_reads_path = os.path.join(path, basename+'_clustered.fa')
        uc_path           = os.path.join(path, basename+'_clustered.uc')
        
        cmd = [
               'usearch', 
               '-cluster_fast', 
               input_reads_path,
               '-id',
               '1',
               '-centroids',
               output_reads_path,
               '-uc',
               uc_path,
               '-sort',
               'size'
               ]
        
        logging.debug('Running command: %s' % ' '.join(cmd))
        subprocess.check_call(cmd)        
        uc_hash = self._read_uc(uc_path)
        
        return output_reads_path, uc_hash
