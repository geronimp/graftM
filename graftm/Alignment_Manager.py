#!/usr/bin/env python
from Bio import SeqIO
from collections import OrderedDict

class Alignment_Manager:
              
    def alignment_correcter(self, alignment_file_list, output_file_name):
        
        corrected_sequences = {}

        for alignment_file in alignment_file_list:
            insert_list = [] # Define list containing inserted positions to be removed (lower case characters)

            sequence_list = list(SeqIO.parse(open(alignment_file, 'r'), 'fasta'))

            for sequence in sequence_list: # For each sequence in the alignment
        
                for idx, nt in enumerate(list(sequence.seq)): # For each nucleotide in the sequence
        
                    if nt.islower(): # Check for lower case character
                        insert_list.append(idx) # Add to the insert list if it is
        
            insert_list = list(OrderedDict.fromkeys(sorted(insert_list, reverse = True))) # Reverse the list and remove duplicate positions
        
            
        
            for sequence in sequence_list: # For each sequence in the alignment
                new_seq = list(sequence.seq) # Define a list of sequences to be iterable list for writing
        
                for position in insert_list: # For each position in the removal list
                    del new_seq[position] # Delete that inserted position in every sequence
                
                corrected_sequences['>'+sequence.id+'\n'] = ''.join(new_seq)+'\n'
        
        
        with open(output_file_name, 'w') as output_file: # Create an open file to write the new sequences to
            for fasta_id, fasta_seq in corrected_sequences.iteritems():
                output_file.write(fasta_id)
                output_file.write(fasta_seq)