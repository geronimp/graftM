#!/usr/bin/env python

import subprocess
import os
import json
import shutil

import graftm.getaxnseq 

from graftm.hmmer import Hmmer
from graftm.messenger import Messenger

class Create:
    
    def __init__(self): 
        self.h=Hmmer(None, None)
        self.the_trash=[]
    
    def buildHmm(self, alignment, base): 
        hmm = base + ".hmm" # Set a name for a hmm
        
        cmd = "hmmbuild --dna %s %s >/dev/null" % (hmm, alignment) # Build the command to build the hmm
        subprocess.check_call(cmd, shell=True) # Call the command
        self.the_trash.append(hmm)
        return hmm
    
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
    
    def buildTree(self, alignment, base): 
        log_file = base + ".tre.log"
        tre_file = base + ".tre"
        
        cmd = "FastTreeMP -quiet -gtr -nt -log %s %s > %s 2>/dev/null" % (log_file, alignment, tre_file)
        subprocess.check_call(cmd, shell=True) # Call the command
        
        self.the_trash += [log_file, tre_file]
        return log_file, tre_file

    def callTaxitCreate(self, base, aln_file, tre, tre_log, tax, seq):
        refpkg = base + ".refpkg"
        
        cmd = "taxit create --quiet -f %s -P %s -t %s -s %s -c -l  %s -T %s -i %s 1>/dev/null" % (aln_file, refpkg, tre, tre_log, base, tax, seq)
        subprocess.check_call(cmd, shell=True)
        
        return refpkg
    
    def compile(self, base, refpkg, hmm, contents): 
        gpkg = self.gpkg_name(base)
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
                
    def gpkg_name(self, base):
        return "%s.gpkg" % base
        
    def main(self, hmm, 
             alignment, 
             sequences, 
             taxonomy, 
             tree, 
             tree_stats,
             sequences_aligned_to_hmm): 
        base=os.path.basename(sequences).split('.')[0]
        
        if os.path.isdir(self.gpkg_name(base)):
            raise Exception("Output directory %s already exists!" % self.gpkg_name(base))    
        
        Messenger().header("Building gpkg for %s" % base)
        # Initially, build the HMM if one is not provided.
        if hmm==None: 
            Messenger().message("Building HMM")
            hmm=self.buildHmm(alignment, base)
        
        # Use the hmm to re-align the sequences if necessary
        if sequences_aligned_to_hmm:
            Messenger().message("Using sequences previously aligned to HMM: %s" % sequences_aligned_to_hmm)
            output_alignment = sequences_aligned_to_hmm
        else:
            Messenger().message("Aligning sequences to HMM")
            output_alignment = self.alignSequences(hmm, sequences, base)
        
        # Build the tree
        if tree and tree_stats:
            log_file = tree_stats
            tre_file = tree
        else:
            Messenger().message("Building tree")
            log_file, tre_file = self.buildTree(output_alignment, base)
        
        # Create tax and seqinfo .csv files
        Messenger().message("Building seqinfo and taxonomy file")
        seq, tax = graftm.getaxnseq.main(base, taxonomy)
        self.the_trash += [seq, tax]
        
        # Create the reference package
        Messenger().message("Creating reference package")
        refpkg = self.callTaxitCreate(base, output_alignment, tre_file, log_file, tax, seq)
        
        # Compile the gpkg
        Messenger().message("Compiling gpkg")
        contents = {"aln_hmm": hmm,
                    "search_hmm": [hmm],
                    "rfpkg": refpkg,
                    "TC":False}
        self.compile(base, refpkg, hmm, contents)

        Messenger().message("Cleaning up")
        self.cleanup(self.the_trash)
        
        Messenger().header("Finished\n")
        return
