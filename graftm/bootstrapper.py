import logging
import tempfile
import shutil
import os
from graftm.hmmer import Hmmer
from graftm.unpack_sequences import UnpackRawReads
import subprocess

class Bootstrapper:
    def generate_hmm_from_contigs(self, search_hmm_files, contig_files, 
                                  maximum_range, threads, evalue,
                                  min_orf_length, restrict_read_length, 
                                  output_hmm_file):
        '''Given a collection of search_hmm_files, search the contigs in 
        contig_files, and generate an HMM from the resulting hits, outputing
        it as output_hmm_file.
        
        Parameters
        ----------
        search_hmm_files: list of str
            list of HMM files to search with
        contig_files: list of str
            list of files to search
        threads, evale, min_orf_length, restrict_read_length:
            as per search_and_extract_orfs_matching_protein_database
        output_hmm_file: str
            path to output file
        
        Returns
        -------
        True if genes were recovered, else False'''
        
        hmmer = Hmmer(search_hmm_files)
        
        with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap_orfs') as orfs:
            logging.info("Finding bootstrap hits in provided contigs..")
            for contig_file in contig_files:
                logging.debug("Finding bootstrap hits in %s.." % contig_file)
                unpack = UnpackRawReads(contig_file)
                
                with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap') as \
                                                        hit_reads_orfs_fasta:
                    # search and extract matching ORFs
                    hmmer.search_and_extract_orfs_matching_protein_database(\
                                                        unpack,
                                                          'hmmsearch',
                                                          maximum_range,
                                                          threads,
                                                          evalue,
                                                          min_orf_length,
                                                          restrict_read_length,
                                                          None,
                                                          None,
                                                          None,
                                                          hit_reads_orfs_fasta.name)
                    # Append to the file
                    shutil.copyfileobj(hit_reads_orfs_fasta.name, orfs.name)
            
            # Now have a fasta file of ORFs.
            # Check to make sure the file is not zero-length
            if os.stat(orfs.name).st_size == 0:
                logging.warn("Failed to find any ORFs in the bootstrap contigs, continuing without bootstrap")
                return False
            
            # Run mafft to align them
            with tempfile.NamedTemporaryFile("graftm_bootstrap_aln") as aln:
                cmd = "mafft %s >%s" % (orfs.name, aln.name)
                logging.info("Aligning bootstrap hits..")
                logging.debug("Running alignment cmd: %s" % cmd)
                subprocess.check_call(cmd, shell=True)
            
                # Run hmmbuild to create an HMM
                cmd = "hmmbuild %s %s" % (output_hmm_file, aln.name)
                logging.info("Building HMM from bootstrap hits..")
                logging.debug("Running cmd: %s" % cmd)
                subprocess.check_call(cmd, shell=True)
                
                return True
            
                