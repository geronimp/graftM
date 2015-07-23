import logging
import tempfile
import shutil
import os
import subprocess

from graftm.hmmer import Hmmer
from graftm.unpack_sequences import UnpackRawReads

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
                    with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap2') as \
                                                        hmmsearch_output_table:
                        with tempfile.NamedTemporaryFile(prefix='graftm_bootstrap3') as \
                                                        hit_reads_fasta:
                            hmmer.search_and_extract_orfs_matching_protein_database(\
                                    unpack,
                                    'hmmsearch',
                                    maximum_range,
                                    threads,
                                    evalue,
                                    min_orf_length,
                                    restrict_read_length,
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
                logging.warn("Failed to find any matching ORFs in the bootstrap contigs, continuing without bootstrap")
                return False
            
            # Run mafft to align them
            with tempfile.NamedTemporaryFile(prefix="graftm_bootstrap_aln") as aln:
                cmd = "mafft --auto %s >%s 2>/dev/null" % (orfs.name, aln.name)
                logging.info("Aligning bootstrap hits..")
                logging.debug("Running alignment cmd: %s" % cmd)
                subprocess.check_call(cmd, shell=True)
            
                # Run hmmbuild to create an HMM
                cmd = "hmmbuild --amino %s %s >/dev/null 2>/dev/null" % (output_hmm_file, aln.name)
                logging.info("Building HMM from bootstrap hits..")
                logging.debug("Running cmd: %s" % cmd)
                subprocess.check_call(cmd, shell=True)
                
                return True
            
                