import logging
import tempfile
import shutil
import os
import extern

from graftm.hmmer import Hmmer
from graftm.unpack_sequences import UnpackRawReads

class Bootstrapper:
    def __init__(self, **kwargs):
        '''
        Parameters
        ----------
        search_hmm_files: list of str
            list of HMM files to search with
        graftm_package: GraftMPackage
            use the search HMMs from this graftm package in addition to
            those specified in search_hmm_files
        threads, evale, min_orf_length, restrict_read_length:
            as per hmmer.search_and_extract_orfs_matching_protein_database
        '''
        self.search_hmm_files = kwargs.pop('search_hmm_files',[])
        self.maximum_range = kwargs.pop('maximum_range',None)
        self.threads = kwargs.pop('threads',None)
        self.evalue = kwargs.pop('evalue',None)
        self.min_orf_length = kwargs.pop('min_orf_length',None)
        graftm_package = kwargs.pop('graftm_package',None)
        if graftm_package:
            for h in graftm_package.search_hmm_files():
                self.search_hmm_files.append(h)
            if self.maximum_range is None:
                self.maximum_range = graftm_package.maximum_range()
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
                    
        
    def generate_hmm_from_contigs(self, contig_files, output_hmm_file):
        '''Given a collection of search_hmm_files, search the contigs in 
        contig_files, and generate an HMM from the resulting hits, outputting
        it as output_hmm_file.
        
        Parameters
        ----------
        contig_files: list of str
            list of files to search
        output_hmm_file: str
            path to output file
        
        Returns
        -------
        True if genes were recovered, else False'''
        
        hmmer = Hmmer(self.search_hmm_files)
        
        with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap_orfs') as orfs:
            logging.info("Finding bootstrap hits in provided contigs..")
            for contig_file in contig_files:
                logging.debug("Finding bootstrap hits in %s.." % contig_file)
                unpack = UnpackRawReads(contig_file)
                
                with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap') as \
                                                        hit_reads_orfs_fasta:
                    # search and extract matching ORFs
                    with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap2') as \
                                                        hmmsearch_output_table:
                        with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap3') as \
                                                        hit_reads_fasta:
                            hmmer.search_and_extract_orfs_matching_protein_database(\
                                    unpack,
                                    'hmmsearch',
                                    self.maximum_range,
                                    self.threads,
                                    self.evalue,
                                    self.min_orf_length,
                                    None,
                                    None,
                                    hmmsearch_output_table.name,
                                    hit_reads_fasta.name,
                                    hit_reads_orfs_fasta.name)
                    # Append to the file
                    shutil.copyfileobj(open(hit_reads_orfs_fasta.name), orfs)
            
            # Now have a fasta file of ORFs.
            # Check to make sure the file is not zero-length
            orfs.flush()
            if os.stat(orfs.name).st_size == 0:
                logging.warn("Failed to find any matching ORFs in the bootstrap contigs")
                return False
            
            # Run mafft to align them
            with tempfile.NamedTemporaryFile(prefix="graftm_bootstrap_aln") as aln:
                cmd = "mafft --auto %s >%s" % (orfs.name, aln.name)
                logging.info("Aligning bootstrap hits..")
                logging.debug("Running alignment cmd: %s" % cmd)
                extern.run(cmd)
            
                # Run hmmbuild to create an HMM
                cmd = "hmmbuild --amino %s %s >/dev/null" % (output_hmm_file, aln.name)
                logging.info("Building HMM from bootstrap hits..")
                logging.debug("Running cmd: %s" % cmd)
                extern.run(cmd)
                
                return True
            
                