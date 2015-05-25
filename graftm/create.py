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
from graftm.housekeeping import HouseKeeping


class Create:
    
    def __init__(self): 
        self.h=Hmmer(None, None)
        self.hk = HouseKeeping()
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
        
        
        cmd = "hmmbuild --dna %s %s >/dev/null" % (hmm, alignment) # Build the command to build the hmm
        
        subprocess.check_call(cmd, shell=True) # Call the command
        return hmm
    
    def pipeType(self, hmm):

        type=[x.split() for x in open(hmm).readlines() if x.startswith('ALPH') or x.startswith('LENG')]
        for item in type:
            if item[0]=='ALPH':       
                if type[1]=='amino':
                    ptype='aa'
                elif type[1]=='DNA' or 'RNA':
                    ptype='na'
                else:
                    raise Exception("Unfamiliar HMM type: %s" % (type[1]))
            elif item[0]=='LENG':
                leng=item[1]
            else:
                raise Exception("Programming Error: Misread HMM file")
        return ptype, leng
    
    def checkAlnLength(self, alignment):
        #seq_format=self.hk.guess_sequence_input_file_format(alignment)
        return len(list(SeqIO.parse(open(alignment, 'r'), 'fasta'))[0].seq)
        
    def alignSequences(self, hmm, sequences, base): 
        stockholm_alignment = base +".aln.sto" # Set an output path for the alignment
        fasta_alignment = base+".insertions.aln.fa" # Set an output path for the alignment
        corrected_fasta_alignment = base+".aln.fa" # Set an output path for the alignment

        cmd = "hmmalign --trim -o %s %s %s" % (stockholm_alignment, hmm, sequences) # Build the command to align the sequences
        subprocess.check_call(cmd, shell=True) # Call the command
        cmd = "seqmagick convert %s %s" % (stockholm_alignment, fasta_alignment)
        subprocess.check_call(cmd, shell=True) # Call the command
        self.h.alignment_correcter([fasta_alignment], corrected_fasta_alignment)
        
        self.the_trash += [stockholm_alignment, corrected_fasta_alignment, fasta_alignment]
        return corrected_fasta_alignment
    
    def buildTree(self, alignment, base, ptype): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        if ptype == 'na': # If it's a nucleotide sequence
            cmd = "FastTreeMP -quiet -gtr -nt -log %s %s > %s 2>/dev/null" % (log_file, alignment, tre_file)
            subprocess.check_call(cmd, shell=True) # Call the command
        else: # Or if its an amino acid sequence
            cmd = "FastTreeMP -quiet -log %s %s > %s 2>/dev/null" % (log_file, alignment, tre_file)
            subprocess.check_call(cmd, shell=True) # Call the command
            
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def callTaxitCreate(self, base, aln_file, tre, tre_log, tax, seq, prefix):
        if prefix:
            refpkg = prefix + ".refpkg"
        else:
            refpkg = base + ".refpkg"
        cmd = "taxit create --quiet -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s 1>/dev/null" % (aln_file, refpkg, tre, tre_log, base, tax, seq)
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
        return

    def cleanup(self, the_trashcan):
        for every_piece_of_junk in the_trashcan:
            if os.path.isdir(every_piece_of_junk):
                shutil.rmtree(every_piece_of_junk)
            else:
                os.remove(every_piece_of_junk)
        
    def main(self, hmm, alignment, sequences, taxonomy, tree, log, prefix): 
        if sequences:
            base=os.path.basename(sequences).split('.')[0]
        else:
            base=os.path.basename(alignment).split('.')[0]
        Messenger().header("Building gpkg for %s" % base)
        # Initially, build the HMM if one is not provided.
        if not tree and not log:
            if hmm and alignment:
                ptype,leng=self.pipeType(hmm)
                if str(self.checkAlnLength(alignment)) != str(leng):
                    logging.info("Alignment length does not match the HMM length, building a new HMM")
                    hmm=self.buildHmm(alignment, base)
                    logging.info("Aligning to HMM built from alignment")
                    output_alignment = self.alignSequences(hmm, sequences, base)
                else:
                    logging.info("Alignment length matches the HMM length, no need to align")
                    output_alignment=alignment
            elif hmm and not alignment:
                logging.info("Using provided HMM to align sequences")
                output_alignment = self.alignSequences(hmm, sequences, base)
            elif alignment and not hmm:
                logging.info("Building HMM from alignment")
                hmm=self.buildHmm(alignment, base)
                logging.info("Aligning to HMM built from alignment")
                output_alignment = self.alignSequences(hmm, sequences, base)
            ptype,leng=self.pipeType(hmm)
            # Build the tree
            logging.info("Building tree")
            log_file, tre_file = self.buildTree(output_alignment, base, ptype)
        else:
            if not tree:
                raise Exception("No Tree provided for log file")
            if not log:
                raise Exception("No log provided for tree")
            try:
                output_alignment=alignment
                tre_file=tree
                log_file=log
            except:
                raise Exception("No alignment file provided")

        # Create tax and seqinfo .csv files
        logging.info("Building seqinfo and taxonomy file")
        seq, tax = graftm.getaxnseq.main(base, taxonomy)
        self.the_trash += [seq, tax]
        
        # Create the reference package
        logging.info("Creating reference package")
        refpkg = self.callTaxitCreate(base, output_alignment, tre_file, log_file, tax, seq, prefix)

        # Compile the gpkg
        logging.info("Compiling gpkg")
        contents = {"aln_hmm": hmm,
                    "search_hmm": [hmm],
                    "rfpkg": refpkg,
                    "TC":False}
        self.compile(base, refpkg, hmm, contents, prefix)

        logging.info("Cleaning up")
        self.cleanup(self.the_trash)
        
        Messenger().header("Finished\n")
        return
