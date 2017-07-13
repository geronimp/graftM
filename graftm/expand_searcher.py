import logging
import tempfile
import shutil
import extern

from graftm.sequence_searcher import SequenceSearcher
from graftm.unpack_sequences import UnpackRawReads
from graftm.sequence_io import SequenceIO

class ExpandSearcher:
    
    # Define constants.
    DIAMOND_SEARCH_METHOD = "diamond"
    HMM_SEARCH_METHOD = "hmmsearch"    

    def __init__(self, **kwargs):
        '''
        Parameters
        ----------
        search_hmm_files: list of str
            list of HMM files to search with
        graftm_package: GraftMPackage
            use the search HMMs from this graftm package in addition to
            those specified in search_hmm_files
        threads, evalue, min_orf_length, restrict_read_length:
            as per hmmer.search_and_extract_orfs_matching_protein_database
        '''
        self.search_hmm_files = kwargs.pop('search_hmm_files',[])
        self.maximum_range = kwargs.pop('maximum_range',None)
        self.threads = kwargs.pop('threads',None)
        self.evalue = kwargs.pop('evalue',None)
        self.min_orf_length = kwargs.pop('min_orf_length',None)
        graftm_package = kwargs.pop('graftm_package',None)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        if graftm_package:
            self.diamond_database = graftm_package.diamond_database_path()
            self.unaligned_sequence_database = graftm_package.unaligned_sequence_database_path()
            if self.search_hmm_files is None:
                self.search_hmm_files = []
            for h in graftm_package.search_hmm_paths():
                if h not in self.search_hmm_files:
                    self.search_hmm_files.append(h)
            if self.maximum_range is None:
                self.maximum_range = graftm_package.maximum_range()
        else:
            self.diamond_database = None
            self.unaligned_sequence_database = None
            
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
                    
        
    def generate_expand_search_database_from_contigs(self, contig_files, output_database_file, 
                                                 search_method):
        '''Given a collection of search_hmm_files, search the contigs in 
        contig_files, and generate an HMM from the resulting hits, outputting
        it as output_database_file.
        
        Parameters
        ----------
        contig_files: list of str
            list of files to search
        output_database_file: str
            path to output file
        search_method: str
            "diamond" or "hmmsearch", to specify search method to use and what
            type of database to build. 
        
        Returns
        -------
        True if genes were recovered, else False'''
        
        ss = SequenceSearcher(self.search_hmm_files)
        seqio = SequenceIO()
        if search_method == self.DIAMOND_SEARCH_METHOD:
            if self.diamond_database == None or self.unaligned_sequence_database == None:
                logging.warning("Cannot expand_search continue with no diamond database or unaligned sequences.") 
                return False
        
        with tempfile.NamedTemporaryFile(prefix='graftm_expand_search_orfs') as orfs:
            logging.info("Finding expand_search hits in provided contigs..")
            for contig_file in contig_files:
                logging.debug("Finding expand_search hits in %s.." % contig_file)
                unpack = UnpackRawReads(contig_file)
                
                with tempfile.NamedTemporaryFile(prefix='graftm_expand_search') as \
                                                        hit_reads_orfs_fasta:
                    # search and extract matching ORFs
                    with tempfile.NamedTemporaryFile(prefix='graftm_expand_search2') as \
                                                        hmmsearch_output_table:
                        with tempfile.NamedTemporaryFile(prefix='graftm_expand_search3') as \
                                                        hit_reads_fasta:
                            ss.search_and_extract_orfs_matching_protein_database(\
                                    unpack,
                                    search_method,
                                    self.maximum_range,
                                    self.threads,
                                    self.evalue,
                                    self.min_orf_length,
                                    None,
                                    (self.diamond_database if self.diamond_database else None),
                                    hmmsearch_output_table.name,
                                    hit_reads_fasta.name,
                                    hit_reads_orfs_fasta.name)
                    # Append to the file
                    shutil.copyfileobj(open(hit_reads_orfs_fasta.name), orfs)
            
            # Now have a fasta file of ORFs.
            # Check to make sure the file is not zero-length
            orfs.flush()
            

            
            with tempfile.NamedTemporaryFile(prefix="graftm_expand_search_aln") as aln:
                    
                if search_method == self.HMM_SEARCH_METHOD:
                    
                    # Check that there is more than one sequence to align.
                    if len(seqio.read_fasta_file(orfs.name)) <= 1:# Just to build on this, you need to check if there is > 1 hit
                                                                  # otherwise mafft will fail to align, causing a crash when hmmbuild is 
                                                                  # run on an empty file.
                        logging.warn("Failed to find two or more matching ORFs in the expand_search contigs")
                        return False
                    
                    # Run mafft to align them
                    cmd = "mafft --auto %s >%s" % (orfs.name, aln.name)
                    logging.info("Aligning expand_search hits..")
                    extern.run(cmd)
                
                    # Run hmmbuild to create an HMM
                    cmd = "hmmbuild --amino %s %s >/dev/null" % (output_database_file, aln.name)
                    logging.info("Building HMM from expand_search hits..")

                    extern.run(cmd)
                
                elif search_method == self.DIAMOND_SEARCH_METHOD:

                    # Concatenate database with existing database
                    with tempfile.NamedTemporaryFile(prefix="concatenated_database") as databasefile:
                        for f in [orfs.name, self.unaligned_sequence_database]:
                            for line in open(f):
                                databasefile.write(line)
                        databasefile.flush()
                                
                        # Run diamond make to create a diamond database
                        cmd = "diamond makedb --in '%s' -d '%s'" % (databasefile.name, output_database_file)
                        logging.info("Building a diamond database from expand_search hits..")
                        extern.run(cmd)
                
                else:
                    raise Exception("Search method not recognised: %s" % search_method)
                    return False
                
                return True
            
                