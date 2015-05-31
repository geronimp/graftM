#!/usr/bin/env python

import subprocess
import os
import json
import shutil
import tempfile
import logging

import graftm.getaxnseq 
from Bio import SeqIO
from graftm.hmmer import Hmmer


class Create:
    
    def __init__(self): 
        self.h=Hmmer(None, None)
        self.the_trash=[]
    
    def buildHmm(self, alignment, base): 
        counter=0
        if os.path.isfile(base + ".hmm"):
            counter=0
            while os.path.isfile(base + ".hmm"):
                base=base+'_%s' % (str(counter))
                counter+=1
            hmm = base + ".hmm"
        else:   
            hmm = base + ".hmm" # Set a name for a hmm
        cmd = "hmmbuild %s %s >/dev/null" % (hmm, alignment) # Build the command to build the hmm
        logging.debug("Calling command: %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        return hmm
    
    def pipeType(self, hmm):
        logging.debug("Setting pipeline type")
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
        logging.debug("Found alignment type as: %s" % leng)
        return ptype, leng
    
    def checkAlnLength(self, alignment):
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)
        
    def alignSequences(self, hmm, sequences, base): 
        stockholm_alignment = base +".aln.sto" # Set an output path for the alignment
        fasta_alignment = base+".insertions.aln.fa" # Set an output path for the alignment
        corrected_fasta_alignment = base+".aln.fa" # Set an output path for the alignment
        cmd = "hmmalign --trim -o %s %s %s" % (stockholm_alignment, hmm, sequences) # Build the command to align the sequences
        logging.debug("Calling command %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        cmd = "seqmagick convert %s %s" % (stockholm_alignment, fasta_alignment)
        logging.debug("Calling command %s" % (cmd))
        subprocess.check_call(cmd, shell=True) # Call the command
        logging.debug("Correcting alignment")
        self.h.alignment_correcter([fasta_alignment], corrected_fasta_alignment)
        self.the_trash += [stockholm_alignment, corrected_fasta_alignment, fasta_alignment]
        logging.debug("Wrote corrected alignments to: %s" % (corrected_fasta_alignment))
        return corrected_fasta_alignment
    
    def buildTree(self, alignment, base, ptype): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == 'na': # If it's a nucleotide sequence
            cmd = "FastTreeMP -quiet -gtr -nt -log %s %s > %s 2>/dev/null" % (log_file, alignment, tre_file)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
        else: # Or if its an amino acid sequence
            cmd = "FastTreeMP -quiet -log %s %s > %s 2>/dev/null" % (log_file, alignment, tre_file)
            logging.debug("Calling command: %s" % (cmd))
            subprocess.check_call(cmd, shell=True) # Call the command
            
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def callTaxitCreate(self, base, aln_file, tre, tre_log, tax, seq, prefix, no_reroot):
        if prefix:
            refpkg = prefix + ".refpkg"
        else:
            refpkg = base + ".refpkg"
        cmd = "taxit create --quiet -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s 1>/dev/null" % (aln_file, refpkg, tre, tre_log, base, tax, seq)
        if no_reroot:
            cmd += ' --no-reroot'
        logging.debug("Calling command: %s" % cmd)
        subprocess.check_call(cmd, shell=True)
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
                
    def get_hmm_and_alignment(self, sequences, alignment, hmm, base):
        '''Return a HMM file and alignment of sequences to that HMM
        
        Returns
        -------
        * HMM file
        * Alignment of sequences to that HMM
        '''
        if hmm and alignment:
            logging.debug("Found HMM and alignment")
            _,leng=self.pipeType(hmm)
            if str(self.checkAlnLength(alignment)) != str(leng):
                logging.info("Alignment length does not match the HMM length, building a new HMM")
                hmm=self.buildHmm(alignment, base)
                logging.info("Aligning to HMM built from alignment")
                output_alignment = self.alignSequences(hmm, sequences, base)
            else:
                logging.info("Alignment length matches the HMM length, no need to align")
                output_alignment=alignment
        elif hmm and not alignment:
            logging.debug("Found HMM but not alignment")
            logging.info("Using provided HMM to align sequences")
            output_alignment = self.alignSequences(hmm, sequences, base)
        elif alignment:
            logging.debug("Building HMM from alignment")
            hmm=self.buildHmm(alignment, base)
            output_alignment = self.alignSequences(hmm, sequences, base)
        else:
            raise Exception("An alignment or HMM is required")
        return hmm, output_alignment
    
    def generate_tree_log_file(self, tree, alignment, output_log_file_path):
        '''Generate the FastTree log file given a tree and the alignment that
        made that tree
        
        Returns
        -------
        Nothing. The log file as parameter is written as the log file.
        '''
        subprocess.check_call(['bash','-c',"FastTree -nome -mllen -intree '%s' -log %s < %s" %\
                                   (tree, output_log_file_path, alignment)])
        
    def main(self, hmm, alignment, sequences, taxonomy, tree, log, prefix, no_reroot): 
        if sequences:
            base=os.path.basename(sequences).split('.')[0]
        else:
            base=os.path.basename(alignment).split('.')[0]
            
        logging.info("Building gpkg for %s" % base)
        
        # Initially, build the HMM if one is not provided.
        hmm, output_alignment = self.get_hmm_and_alignment(sequences, alignment, hmm, base)
        
        if tree:
            tre_file=tree
            logging.debug("Found tree file")
            if log:
                # User specified a log file, go with that
                log.debug("Using user-specified log file %s" % log)
                log_file = log
            else:
                log.debug("Generating log file")
                log_file = tempfile.NamedTemporaryFile(suffix='.log', prefix='graftm') 
                self.generate_tree_log_file(tree, alignment,
                                            log_file)
        else:
            logging.debug("No tree provided")
            logging.info("Building tree")
            ptype,_ = self.pipeType(hmm)
            log_file, tre_file = self.buildTree(output_alignment, base, ptype)
            
        # Create tax and seqinfo .csv files
        logging.info("Building seqinfo and taxonomy file")
        seq, tax = graftm.getaxnseq.main(base, taxonomy)
        self.the_trash += [seq, tax]
        
        # Create the reference package
        logging.info("Creating reference package")
        refpkg = self.callTaxitCreate(base, output_alignment, tre_file, log_file, tax, seq, prefix, no_reroot)

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
        return
