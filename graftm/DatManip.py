#!/usr/bin/env python
from graftm.Messenger import Messenger

import subprocess
import re
import json
# Constants
FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class DatManip:
    
    def jplace_split(self, jplace_file, alias_hash):
        placement_file = json.load(open(jplace_file))
        jplace_path_list = []
        
        for placement in placement_file['placements']:
            
            hash = {}
            for alias in alias_hash:
                hash = {'p': placement['p'],
                        'nm': [nm for nm in placement['nm'] if nm[0].split('_')[-1] == alias]}
                alias_hash[alias]['place'].append(hash)
         
                 
        # Parse placements
        for alias in alias_hash:

            output = {'fields': placement_file['fields'],
                      'version': placement_file['version'],
                      'tree':  placement_file['tree'],
                      'placements': alias_hash[alias]['place'],
                      'metadata': placement_file['metadata']}
            
            base = alias_hash[alias]['filename'] + '/' + alias_hash[alias]['filename'] + '_placements.jplace'
            
            with open(base, 'w') as out:
                json.dump(output, out, ensure_ascii=False)
            jplace_path_list.append(base)
        
        return jplace_path_list

    
    def csv_to_titles(self, hmm_table_list, graftm_pipeline, readnames_output_path, base):
        '''process hmmsearch/nhmmer results into a list of matching reads/ORFs for D/P respectively, to *_readnames.txt
        '''
        rev_true = False
        evals = {}
        titles_list = []
        reads_list = []
        write_list = []
        title_count = 0
        orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
    
        for hmm_table in hmm_table_list:
    
            for line in open(hmm_table):
    
                if line.startswith('#'):
                    continue
                
                title_count += 1
                
                read_name = str(line.split()[0])
                e_value = line.split()[12]
                direction = line.split()[11]
                
                if direction == '-':
                    rev_true = True
                
                if graftm_pipeline == 'D':
    
                    titles_list.append(read_name)
                    evals[read_name] = [e_value, direction]
    
    
                elif graftm_pipeline == 'P':
                    # The original reads file contains sequences like
                    # >eg and comment
                    # where orfm gives orfs in the form of
                    # >eg_1_2_3 and comment
                    # The read_name here is the orfm style, we want to add to the titles_list
                    # the original form
                    regex_match = orfm_regex.match(read_name)
    
                    if regex_match is not None:
                        titles_list.append(regex_match.groups(0)[0])
    
                    else:
                        raise Exception("Unexpected form of ORF name found: %s" % read_name)
    
                else:
                    raise Exception("Programming error")
            
            reads_list.append(titles_list)
            titles_list = []
        
        # Check if there arereads, and exit if non are found
        if title_count == 0:  
            Messenger().message('0 Reads found! Exiting')
            exit(0)  # and exit
    
    
             
        
        if len(reads_list) == 2:
            for_r = set(reads_list[0])
    
            for read in reads_list[1]:
    
                if read in for_r:
                    write_list.append(read)
    
        elif len(reads_list) == 1:
            write_list = reads_list[0]
    
        Messenger().message('Found %s read(s) in %s' % (len(write_list), base))
        
        with open(readnames_output_path, 'w') as output_file:
    
            for name in write_list:
        
                output_file.write(name + '\n')
        
    
        return evals, len(write_list), rev_true
    
    def extract_from_raw_reads(self, raw_seq_file, hmm, name_file, input_file_format, outdir, sequence_file_name, raw_orf_title, hmm_out_title, orf_titles_path, orf_hmm_out_title, dna_pipe):
        # Run fxtract to obtain reads form original sequence file
        fxtract_cmd = "fxtract -H -X -f %s " % name_file
        if input_file_format == FORMAT_FASTA:
            cmd = fxtract_cmd + " " + raw_seq_file + " > " + sequence_file_name + ' '
            # log
            subprocess.check_call(cmd, shell=True)
        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = fxtract_cmd + raw_seq_file + " | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > " + sequence_file_name + " "
            # log
            subprocess.check_call(fxtract_cmd + raw_seq_file + " | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > " + sequence_file_name + " ", shell=True)
        else:
            raise Exception("Programming error")
    
        # Exit if in the dna pipeline
        if dna_pipe:
            return
    
        # Call orfs on the sequences
        cmd = 'orfm ' + sequence_file_name + ' > ' + raw_orf_title
        # log
        subprocess.check_call(cmd, shell=True)
    
        # Search for the correct reading fram
        cmd = "hmmsearch -o /dev/null --tblout " + hmm_out_title + " " + hmm + " " + raw_orf_title
        # log
        subprocess.check_call(cmd, shell=True)
    
        raw_titles = []
    
        for line in open(hmm_out_title):
    
            if line.startswith('#'):
                continue
    
            split = line.split(' ', 1)
            raw_titles.append(split[0])
    
         
    
        title_file_open = open(orf_titles_path, 'w')
    
        for title in raw_titles:
            title_file_open.write(str(title) + '\n')
    
        title_file_open.close()
        cmd = 'fxtract -H -X -f ' + orf_titles_path + ' ' + raw_orf_title + ' > ' + orf_hmm_out_title
        # log
        subprocess.check_call(cmd, shell=True)
     