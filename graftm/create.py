#!/usr/bin/env python

import subprocess
import os
import json
import shutil
import tempfile
import logging
import subprocess32

from Bio import SeqIO
from graftm.hmmer import Hmmer
from graftm.tree_cleaner import TreeCleaner
from graftm.taxonomy_extractor import TaxonomyExtractor
from graftm.getaxnseq import Getaxnseq
from skbio.tree import TreeNode

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
        
    def get_hmm_and_alignment(self, alignment, base):
        '''Return a HMM file and alignment of sequences to that HMM
        
        Returns
        -------
        * HMM file
        * Alignment of sequences to that HMM
        '''
        logging.debug("Building HMM from alignment")
        hmm, output_alignment = self.buildHmm(alignment, base)
        self.the_trash += [hmm]
        return hmm, output_alignment
    
    def buildHmm(self, alignment, base):
        '''Make a HMM called base.hmm with the given alignment file
        
        Returns
        -------
        A list of two:
        1. Path to the created HMM
        2. Alignment of sequences to this HMM'''
        counter=0
        if os.path.isfile(base + ".hmm"):
            counter=0
            while os.path.isfile(base + ".hmm"):
                base=base+'_%s' % (str(counter))
                counter+=1
            hmm = base + ".hmm"
        else:   
            hmm = base + ".hmm" # Set a name for a hmm

        sto = tempfile.NamedTemporaryFile(suffix='.sto',prefix='graftm')
        tempaln = tempfile.NamedTemporaryFile(suffix='.fasta',prefix='graftm')
        cmd = "hmmbuild -O %s '%s' '%s'" % (sto.name,
                                              hmm,
                                              alignment)
        logging.debug("Calling command: %s" % cmd)
        out = subprocess.check_output(cmd, shell=True) # Call the command
        logging.debug("Got STDOUT from hmmbuild: %s" % out)

        cmd = "seqmagick convert --input-format stockholm %s %s" % (sto.name,
                                                              tempaln.name)
        logging.debug("Calling command: %s" % cmd)
        out = subprocess.check_output(['bash','-c',cmd]) # Call the command
        logging.debug("Got STDOUT from seqmagick: %s" % out)
        
        alignment_path = "%s.aln.fa" % base
        self.h.alignment_correcter([tempaln.name], alignment_path)

        return hmm, alignment_path
    
    def pipeType(self, hmm):
        logging.debug("Setting pipeline type..")
        type=[x.split() for x in open(hmm).readlines() if x.startswith('ALPH') or x.startswith('LENG')]
        for item in type:
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
    
    def checkAlnLength(self, alignment):
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)
    
    def buildTree(self, alignment, base, ptype): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == 'na': # If it's a nucleotide sequence
            cmd = "FastTreeMP -quiet -gtr -nt -log %s -out %s %s" % (log_file, tre_file, alignment)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
        else: # Or if its an amino acid sequence
            cmd = "FastTreeMP -quiet -log %s -out %s %s" % (log_file, tre_file, alignment)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
            
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def callTaxitCreate(self, base, aln_file, tre, tre_log, tax, seq, prefix, no_reroot):
        if prefix:
            refpkg = prefix + ".refpkg"
        else:
            refpkg = base + ".refpkg"
            
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
                logging.error('''taxit create failed to run in a small amount of time suggesting that
rerooting was unsuccessful. Unfortunately this tree will need to be rerooted 
manually yourself using a tree editor such as ARB or FigTree.
Once you have a rerooted newick format tree, rerun graftm create
specifying the new tree with --rerooted_tree. The tree file to be rerooted is \'%s\'

''' % tre)
                exit(2)
        return refpkg
    
    def compile(self, base, refpkg, hmm, contents, prefix): 
        if prefix:
            gpkg = prefix + ".gpkg"
        else:
            gpkg = base + ".gpkg"
        if os.path.isdir(gpkg): 
            raise Exception("Detected gpkg with name %s" % (gpkg))
        os.mkdir(gpkg)
        os.rename(refpkg, os.path.join(gpkg, refpkg))
        shutil.copyfile(hmm, os.path.join(gpkg, os.path.basename(hmm)))
        json.dump(contents, open(os.path.join(gpkg, 'CONTENTS.json'), 'w'))

    def cleanup(self, the_trashcan):
        for every_piece_of_junk in the_trashcan:
            if os.path.isdir(every_piece_of_junk):
                shutil.rmtree(every_piece_of_junk)
            else:
                os.remove(every_piece_of_junk)
    
    def generate_tree_log_file(self, tree, alignment, output_tree_file_path,
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
        logging.debug("Running log creation command %s" % cmd)
        subprocess.check_call(['bash','-c',cmd])
        
    def main(self, alignment, **kwargs):
        taxonomy = kwargs.pop('taxonomy',None)
        rerooted_tree = kwargs.pop('rerooted_tree',None)
        tree_log = kwargs.pop('tree_log', None)
        prefix = kwargs.pop('prefix', None)
        rerooted_annotated_tree = kwargs.pop('rerooted_annotated_tree', None)
        min_aligned_percent = kwargs.pop('min_aligned_percent',0.0)
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        base=os.path.basename(alignment).split('.')[0]
            
        logging.info("Building gpkg for %s" % base)
        
        # align sequences to HMM (and potentially build hmm from alignment)
        hmm, output_alignment = self.get_hmm_and_alignment(alignment, base)
        ptype,hmm_length = self.pipeType(hmm)
        import IPython ; IPython.embed()
        logging.info("Checking for incorrect or fragmented reads")
        insufficiently_aligned_sequences = self._check_reads_hit(open(output_alignment),
                                                                 min_aligned_percent)
        if len(insufficiently_aligned_sequences) > 0:
            logging.error("One or more alignments do not span > %.2f %% of HMM" % (min_aligned_percent*100))
            for s in insufficiently_aligned_sequences:
                logging.error("Insufficient alignment of %s" % s)
            raise Exception("One or more sequences did not span a sufficient"+
                " amount of the alignment suggesting that potentially they should "
                "not be included in the alignment e.g. '%s'" % insufficiently_aligned_sequences[0])
        logging.debug("Found no sequences of insufficient length")
            
        # Create tree unless one was provided
        if not rerooted_tree and not rerooted_annotated_tree:
            logging.debug("No tree provided")
            logging.info("Building tree")
            log_file, tre_file = self.buildTree(output_alignment, base, ptype)
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
            TreeCleaner().match_alignment_and_tree_sequence_ids(output_alignment,
                                                                tre_file)
            if tree_log:
                # User specified a log file, go with that
                logging.debug("Using user-specified log file %s" % tree_log)
                log_file = tree_log
            else:
                logging.info("Generating log file")
                log_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree_log', prefix='graftm') 
                log_file = log_file_tempfile.name
                input_tree_file = tre_file
                tre_file1_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file1 = tre_file1_tempfile.name
                # Make the newick file simple (ie. un-arb it) for fasttree
                TreeCleaner().clean_newick_file_for_fasttree_input(input_tree_file, tre_file1)
                tre_file_tempfile = tempfile.NamedTemporaryFile(suffix='.tree', prefix='graftm')
                tre_file = tre_file_tempfile.name
                self.generate_tree_log_file(tre_file1, output_alignment,
                                            tre_file, log_file, ptype)
            
        # Create tax and seqinfo .csv files
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
            raise Exception("Taxonomy is required somehow e.g. by --taxonomy or --rerooted_annotated_tree")
        gtns.write_taxonomy_and_seqinfo_files(taxonomy_definition,
                                              tax,
                                              seq)
        
        # Create the reference package
        logging.info("Creating reference package")
        refpkg = self.callTaxitCreate(base, output_alignment, tre_file, 
                                      log_file, tax, seq, prefix, no_reroot)

        # Compile the gpkg
        logging.info("Compiling gpkg")
        contents = {"aln_hmm": hmm,
                    "search_hmm": [hmm],
                    "rfpkg": refpkg,
                    "TC":False}
        self.compile(base, refpkg, hmm, contents, prefix)

        logging.info("Cleaning up")
        self.cleanup(self.the_trash)
        
        logging.info("Finished\n")
