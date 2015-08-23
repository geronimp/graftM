#!/usr/bin/env python

import os
import shutil
import tempfile
import logging
import subprocess32
import extern

from Bio import SeqIO
from graftm.hmmer import Hmmer
from graftm.tree_cleaner import TreeCleaner
from graftm.taxonomy_extractor import TaxonomyExtractor
from graftm.getaxnseq import Getaxnseq
from graftm.deduplicator import Deduplicator
from skbio.tree import TreeNode
from graftm.sequence_io import SequenceIO
import tempdir
from graftm.graftm_package import GraftMPackageVersion2

class Create:
    
    def __init__(self): 
        self.h=Hmmer(None, None)
        self.the_trash=[]
    
    def _check_reads_hit(self, alignment_io, min_aligned_fraction):
        '''Given an alignment return a list of sequence names that are less
        than the min_aligned_fraction'''
        to_return = []
        alignment_length = None
        for s in SeqIO.parse(alignment_io, "fasta"):
            if not alignment_length:
                alignment_length = len(s.seq)
                min_length = int(min_aligned_fraction * alignment_length)
                logging.debug("Determined min number of aligned bases to be %s" % min_length)
            elif len(s.seq) != alignment_length:
                raise Exception("Alignment file appears to not be of uniform length")

            num_unaligned = s.seq.count('-')
            num_aligned = alignment_length-num_unaligned
            logging.debug("Sequence %s has %d aligned positions" % (s.name, alignment_length-num_unaligned))
            if num_aligned <= min_length:
                to_return.append(s.name)
        return to_return
        
    def _get_hmm_from_alignment(self, alignment, hmm_filename, output_alignment_filename):
        '''Return a HMM file and alignment of sequences to that HMM
        
        Parameters
        ----------
        alignment: str
            path to aligned proteins
        hmm_filename: str
            write the hmm to this file path
        output_alignment_filename: str
            write the output alignment to this file path
        
        Returns
        -------
        Nothing
        '''
        logging.debug("Building HMM from alignment")
        
        sto = tempfile.NamedTemporaryFile(suffix='.sto',prefix='graftm')
        tempaln = tempfile.NamedTemporaryFile(suffix='.fasta',prefix='graftm')
        cmd = "hmmbuild -O %s '%s' '%s'" % (sto.name,
                                              hmm_filename,
                                              alignment)
        out = extern.run(cmd)
        logging.debug("Got STDOUT from hmmbuild: %s" % out)

        cmd = "seqmagick convert --input-format stockholm %s %s" % (sto.name,
                                                              tempaln.name)
        out = extern.run(cmd)
        logging.debug("Got STDOUT from seqmagick: %s" % out)
        
        Hmmer(hmm_filename).alignment_correcter([tempaln.name], output_alignment_filename)
        
    def _align_sequences_to_hmm(self, hmm_file, sequences_file, output_alignment_file):
        '''Align sequences to an HMM, and return a path to an alignment of
        these proteins after cleanup so that they can be used for tree-making
        
        Parameters
        ----------
        sequences_file: str
            path to file of unaligned protein sequences
        hmm_file: str
            path to hmm file
        output_alignment_file: str
            write alignment to this file
        
        Returns
        -------
        nothing
        '''
        hmmer = Hmmer(hmm_file)
        tempalign = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.fasta')
        hmmer.hmmalign_sequences(hmm_file, sequences_file, tempalign.name)
        hmmer.alignment_correcter([tempalign.name], output_alignment_file)
    
    def _pipe_type(self, hmm):
        logging.debug("Setting pipeline type..")
        hmm_type=[x.split() for x in open(hmm).readlines() if x.startswith('ALPH') or x.startswith('LENG')]
        for item in hmm_type:
            if item[0]=='ALPH':
                if item[1]=='amino':
                    ptype='aa'
                elif item[1]=='DNA' or 'RNA':
                    ptype='na'
                else:
                    raise Exception("Unfamiliar HMM type: %s" % (item[1]))
            elif item[0]=='LENG':
                leng=item[1]
            else:
                raise Exception("Programming Error: Misread HMM file")
        logging.debug("Set pipeline type as: %s " % ptype)
        logging.debug("Found alignment length as: %s" % leng)
        return ptype, leng
    
    def _check_aln_length(self, alignment):
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)
    
    def _build_tree(self, alignment, base, ptype): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == 'na': # If it's a nucleotide sequence
            cmd = "FastTreeMP -quiet -gtr -nt -log %s -out %s %s" % (log_file, tre_file, alignment)
            extern.run(cmd)
        else: # Or if its an amino acid sequence
            cmd = "FastTreeMP -quiet -log %s -out %s %s" % (log_file, tre_file, alignment)
            extern.run(cmd)
            
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def _taxit_create(self, base, aln_file, tre, tre_log, tax, seq, refpkg, no_reroot):
        cmd = "taxit create -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s"\
            % (aln_file, refpkg, tre, tre_log, base, tax, seq)
            
        if no_reroot:
            cmd += ' --no-reroot'
            logging.debug("Calling command assuming pre-rerooting: %s" % cmd)
            subprocess32.check_call(cmd, shell=True)
        else:
            logging.debug("Calling command: %s" % cmd)
            logging.info("Attempting to run taxit create with rerooting capabilities")
            try:
                subprocess32.check_call(cmd, shell=True, timeout=20)
            except (subprocess32.TimeoutExpired, subprocess32.CalledProcessError):
                nontemp_tax_file = 'graftm_create_taxonomy.csv'
                shutil.copy(tax, nontemp_tax_file)
                nontemp_aln_file = 'graftm_create_alignment.faa'
                shutil.copy(aln_file, nontemp_aln_file)
                
                logging.error('''
taxit create failed to run in a small amount of time suggesting that
rerooting was unsuccessful. Unfortunately this tree will need to be rerooted 
manually yourself using a tree editor such as ARB or FigTree.
Once you have a rerooted newick format tree, rerun graftm create
specifying the new tree with --rerooted_tree. The tree file to be rerooted is \'%s\'''' % tre)

                logging.error('''
When rerunning, please use the following flags for the command line to account
for the fact that some sequences may have been removed during the deduplication
process.

graftM create --taxonomy '%s' --alignment '%s' aln_file

(plus other relevant arguments).
''' % (nontemp_tax_file, nontemp_aln_file))
                exit(2)
        return refpkg
    


    def _cleanup(self, the_trashcan):
        for every_piece_of_junk in the_trashcan:
            if os.path.isdir(every_piece_of_junk):
                shutil.rmtree(every_piece_of_junk)
            else:
                os.remove(every_piece_of_junk)
    
    def _define_range(self, sequences):
        '''
        define_range - define the maximum range within which two hits in a db 
        search can be linked. This is defined as 1.5X the average length of all
        reads in the database. 
        
        Parameters
        ----------
        sequences : str
            A path to the sequences in FASTA format. This fasta file is assumed 
            to be in the correct format. i.e. headers start with '>'
            
        Returns
        -------
        max_range : int 
            As described above, 1.5X the size of the average length of genes 
            within the database
        '''
        sequence_count = 0
        total_sequence = 0

        for record in SeqIO.parse(open(sequences), 'fasta'): # For every sequence in the output 
            total_sequence+=1 # increment the total sequence count
            sequence_count+=len(record.seq) # add the length of that sequence to the total string count
                
        # Get the average, multiply by 1.5, and return
        max_range = (sequence_count/total_sequence)*1.5
        return max_range 

    def _generate_tree_log_file(self, tree, alignment, output_tree_file_path,
                               output_log_file_path, residue_type):
        '''Generate the FastTree log file given a tree and the alignment that
        made that tree
        
        Returns
        -------
        Nothing. The log file as parameter is written as the log file.
        '''
        if residue_type=='na':
            cmd = "FastTree -quiet -gtr -nt -nome -mllen -intree '%s' -log %s -out %s %s" %\
                                       (tree, output_log_file_path,
                                        output_tree_file_path, alignment)
        elif residue_type=='aa':
            cmd = "FastTree -quiet -nome -mllen -intree '%s' -log %s -out %s %s" %\
                                       (tree, output_log_file_path, 
                                        output_tree_file_path, alignment)
        extern.run(cmd)
    
    def _remove_sequences_from_alignment(self, sequence_names, input_alignment_file, output_alignment_file):
        '''Remove sequences from the alignment file that have names in
        sequence_names
        
        Parameters
        ----------
        sequence_names: list of str
            names of sequences to remove
        input_alignment_file: str
            path to alignment file to remove from
        output_alignment_file: str
            path to alignment file to write to
            
        Returns
        -------
        int: number of sequences written to file'''
        nameset = set(sequence_names)
        num_written = 0
        with open(output_alignment_file, 'w') as f:
            for s in SeqIO.parse(open(input_alignment_file), "fasta"):
                if s.name not in nameset:
                    f.write(">"+s.name+"\n")
                    f.write(str(s.seq)+"\n")
                    num_written += 1
        return num_written
                

        
    def main(self, **kwargs):
        alignment = kwargs.pop('alignment',None)
        sequences = kwargs.pop('sequences',None)
        taxonomy = kwargs.pop('taxonomy',None)
        rerooted_tree = kwargs.pop('rerooted_tree',None)
        tree_log = kwargs.pop('tree_log', None)
        prefix = kwargs.pop('prefix', None)
        rerooted_annotated_tree = kwargs.pop('rerooted_annotated_tree', None)
        user_hmm = kwargs.pop('hmm', None)
        min_aligned_percent = kwargs.pop('min_aligned_percent',0.01)
        taxtastic_taxonomy = kwargs.pop('taxtastic_taxonomy', None)
        taxtastic_seqinfo = kwargs.pop('taxtastic_seqinfo', None)
        force_overwrite = kwargs.pop('force',False)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        locus_name = os.path.basename(sequences).split('.')[0]
        tmp = tempdir.TempDir()
        base = os.path.join(tmp.name, locus_name)
            
        if prefix:
            output_gpkg_path = prefix
        else:
            output_gpkg_path = "%s.gpkg"
            
        if os.path.exists(output_gpkg_path):
            if force_overwrite:
                logging.warn("Deleting previous directory %s" % output_gpkg_path)
                shutil.rmtree(output_gpkg_path)
            else:
                raise Exception("Cowardly refusing to overwrite gpkg to already existing %s" % output_gpkg_path)
        logging.info("Building gpkg for %s" % output_gpkg_path)
        
        # align sequences to HMM (and potentially build hmm from alignment)
        output_alignment_f = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.faa')
        output_alignment = output_alignment_f.name
        if user_hmm:
            hmm = user_hmm
            self._align_sequences_to_hmm(hmm, sequences, output_alignment)
        else:
            hmm_f = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.hmm')
            hmm = hmm_f.name
            self._get_hmm_from_alignment(alignment, hmm, output_alignment)
        ptype, _ = self._pipe_type(hmm)
        logging.info("Checking for incorrect or fragmented reads")
        insufficiently_aligned_sequences = self._check_reads_hit(open(output_alignment),
                                                                 min_aligned_percent)
        if len(insufficiently_aligned_sequences) > 0:
            logging.warn("One or more alignments do not span > %.2f %% of HMM" % (min_aligned_percent*100))
            for s in insufficiently_aligned_sequences:
                logging.warn("Insufficient alignment of %s, not including this sequence" % s)
            output_alignment_f2 = tempfile.NamedTemporaryFile(prefix='graftm', suffix='.aln.faa')
            output_alignment2 = output_alignment_f2.name
            num_sequences = self._remove_sequences_from_alignment(insufficiently_aligned_sequences, output_alignment,output_alignment2)
            output_alignment = output_alignment2
            logging.info("After removing insufficiently aligned sequences, left with %i sequences aligned" % num_sequences)
            if num_sequences < 4:
                raise Exception("Too few sequences remaining in alignment after removing insufficiently aligned sequences: %i" % num_sequences)
        else:
            logging.debug("Found no sequences of insufficient length")
        
        # Read in taxonomy somehow
        gtns = Getaxnseq()
        if rerooted_annotated_tree:
            logging.info("Building seqinfo and taxonomy file from input annotated tree")
            taxonomy_definition = TaxonomyExtractor().taxonomy_from_annotated_tree(\
                    TreeNode.read(open(rerooted_annotated_tree)))
        elif taxonomy:
            logging.info("Building seqinfo and taxonomy file from input taxonomy")
            taxonomy_definition = gtns.read_taxonomy_file(taxonomy)
        elif taxtastic_seqinfo and taxtastic_taxonomy:
            logging.info("Reading taxonomy from taxtastic taxonomy and seqinfo files")
            taxonomy_definition = gtns.read_taxtastic_taxonomy_and_seqinfo\
                (open(taxtastic_taxonomy), 
                 open(taxtastic_seqinfo))
        else:
            raise Exception("Taxonomy is required somehow e.g. by --taxonomy or --rerooted_annotated_tree")
        
        # Make sure each sequence has been assigned a taxonomy:
        seqio = SequenceIO()
        aligned_sequence_objects = seqio.read_fasta_file(output_alignment)
        unannotated = []
        for s in aligned_sequence_objects:
            if s.name not in taxonomy_definition:
                unannotated.append(s.name)
        if len(unannotated) > 0:
            for s in unannotated:
                logging.error("Unable to find sequence '%s' in the taxonomy definition" % s)
            raise Exception("All sequences must be assigned a taxonomy, cannot continue")
        
        # Deduplicate sequences - pplacer cannot handle these
        logging.info("Deduplicating sequences")
        dedup = Deduplicator()
        deduplicated_arrays = dedup.deduplicate(aligned_sequence_objects)
        deduplicated_taxonomy = dedup.lca_taxonomy(deduplicated_arrays, taxonomy_definition)
        deduplicated_taxonomy_hash = {}
        for i, tax in enumerate(deduplicated_taxonomy):
            deduplicated_taxonomy_hash[deduplicated_arrays[i][0].name] = tax
        deduplicated_alignment_file = base+"_deduplicated_aligned.fasta"
        seqio.write_fasta_file([seqs[0] for seqs in deduplicated_arrays],
                               deduplicated_alignment_file)
        logging.info("Removed %i sequences as duplicates, leaving %i non-identical sequences"\
                     % ((len(aligned_sequence_objects)-len(deduplicated_arrays)),
                        len(deduplicated_arrays)))
            
        # Create tree unless one was provided
        if not rerooted_tree and not rerooted_annotated_tree:
            logging.debug("No tree provided")
            logging.info("Building tree")
            log_file, tre_file = self._build_tree(deduplicated_alignment_file, base, ptype)
            no_reroot = False
        else:
            if rerooted_tree:
                logging.debug("Found unannotated pre-rerooted tree file %s" % rerooted_tree)
                tre_file=rerooted_tree
            elif rerooted_annotated_tree:
                logging.debug("Found annotated pre-rerooted tree file %s" % rerooted_tree)
                tre_file=rerooted_annotated_tree
            else:
                raise
            no_reroot = True
            
            # Remove any sequences from the tree that are duplicates
            cleaner = TreeCleaner()
            tree = TreeNode.read(open(tre_file))
            removed_sequence_names = []
            for group in deduplicated_arrays:
                [removed_sequence_names.append(s.name) for s in group[1:]]        
            cleaner.remove_sequences(tree, removed_sequence_names)
            
            # Ensure there is nothing amiss now as a user-interface thing
            cleaner.match_alignment_and_tree_sequence_ids(\
                [g[0].name for g in deduplicated_arrays], tree)
            
            if tree_log:
                # User specified a log file, go with that
                logging.debug("Using user-specified log file %s" % tree_log)
                log_file = tree_log
            else:
                logging.info("Generating log file")
                log_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree_log', prefix='graftm') 
                log_file = log_file_tempfile.name
                # Make the newick file simple (ie. un-arb it) for fasttree
                cleaner.clean_newick_file_for_fasttree_input(tree)
                tre_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file = tre_file_tempfile.name
                with tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm') as f:
                    tree.write(f)
                    f.flush()
                    self._generate_tree_log_file(f.name, deduplicated_alignment_file,
                                            tre_file, log_file, ptype)
            
        # Create tax and seqinfo .csv files
        refpkg = "%s.refpkg" % output_gpkg_path
        self.the_trash.append(refpkg)
        if taxtastic_taxonomy and taxtastic_seqinfo:
            logging.info("Creating reference package")
            refpkg = self._taxit_create(base, deduplicated_alignment_file, tre_file, 
                                          log_file, taxtastic_taxonomy,
                                          taxtastic_seqinfo, refpkg, no_reroot)
        else:
            gtns = Getaxnseq()
            seq = base+"_seqinfo.csv"
            tax = base+"_taxonomy.csv"
            self.the_trash += [seq, tax]
            if rerooted_annotated_tree:
                logging.info("Building seqinfo and taxonomy file from input annotated tree")
                taxonomy_definition = TaxonomyExtractor().taxonomy_from_annotated_tree(\
                        TreeNode.read(open(rerooted_annotated_tree)))
            elif taxonomy:
                logging.info("Building seqinfo and taxonomy file from input taxonomy")
                taxonomy_definition = gtns.read_taxonomy_file(taxonomy)
                
            else:
                raise Exception("Programming error: Taxonomy is required somehow e.g. by --taxonomy or --rerooted_annotated_tree")
            gtns.write_taxonomy_and_seqinfo_files(taxonomy_definition,
                                                  tax,
                                                  seq)
            
            # Create the reference package
            logging.info("Creating reference package")
            refpkg = self._taxit_create(base, deduplicated_alignment_file, tre_file, 
                                          log_file, tax, seq, refpkg, no_reroot)
        
        # Run diamond makedb
        logging.info("Creating diamond database")
        cmd = "diamond makedb --in '%s' -d '%s'" % (sequences, base)
        extern.run(cmd)
        diamondb = '%s.dmnd' % base

        # Get range
        max_range=self._define_range(sequences)
        
        # Compile the gpkg
        logging.info("Compiling gpkg")
        GraftMPackageVersion2.compile(output_gpkg_path, locus_name, refpkg, hmm, diamondb, max_range)

        logging.info("Cleaning up")
        self._cleanup(self.the_trash)
        
        logging.info("Finished\n")
