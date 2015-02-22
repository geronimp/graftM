#!/usr/bin/env python
import subprocess
from graftm.GraftMFiles import GraftMFiles
from graftm.Messenger import Messenger
from Bio import SeqIO


FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class Hmmer:
    
    def __init__(self, hmm):
        self.hmm = hmm
        
    def hmmalign(self, base, sequencefile, summary_dict):

        GMF = GraftMFiles(base)
        for_file = GMF.output_for_path(base)
        rev_file = GMF.output_rev_path(base)
        for_sto_file = GMF.sto_for_output_path(base)
        rev_sto_file = GMF.sto_rev_output_path(base)
        for_conv_file = GMF.conv_output_for_path(base)
        rev_conv_file = GMF.conv_output_rev_path(base)
        ## TODO look at old graftM and re-paste code. for HMMalign. I screwed somthing up in here that makes it twice as slow....
        # If there are reverse complement reads
        if summary_dict['rev_true']:

            evals = summary_dict['evals']

            reverse = []
            forward = []
            
            records = list(SeqIO.parse(open(sequencefile), 'fasta'))
            
            # Split the reads into reverse and forward lists
            for record in records:
                    
                if evals[record.id][1] == '+':
                    forward.append(record)
                    
                elif evals[record.id][1] == '-':
                    reverse.append(record)
                    
                else:
                    Messenger().error_message('Programming error: hmmalign')
                    exit(1)
                
                
            # Write reverse complement and forward reads to files
            with open(for_file, 'w') as for_aln:
                for record in forward:
                    for_aln.write('>'+record.id+'\n')
                    for_aln.write(str(record.seq)+'\n')
                
            with open(rev_file, 'w') as rev_aln:
                for record in reverse:
                    rev_aln.write('>'+record.id+'\n')
                    rev_aln.write(str(record.seq.reverse_complement())+'\n')
            
                            
    
            # HMMalign and convert to fasta format
            cmd = 'hmmalign --trim -o %s %s %s 2>/dev/null; seqmagick convert %s %s' % (for_sto_file, 
                                                                                        self.hmm, 
                                                                                        for_file, 
                                                                                        for_sto_file, 
                                                                                        for_conv_file)
                # log
            subprocess.check_call(cmd, shell=True)
                    
            cmd = 'hmmalign --trim -o %s %s %s 2>/dev/null; seqmagick convert %s %s' % (rev_sto_file, 
                                                                                        self.hmm, 
                                                                                        rev_file, 
                                                                                        rev_sto_file, 
                                                                                        rev_conv_file)
                # log
            subprocess.check_call(cmd, shell=True)
                
            # If there are only forward reads, just hmmalign and be done with it.
        else:
                
            cmd = 'hmmalign --trim -o %s %s %s ; seqmagick convert %s %s' % (for_sto_file, 
                                                                             self.hmm, 
                                                                             sequencefile, 
                                                                             for_sto_file, 
                                                                             for_conv_file)
                # log
            subprocess.check_call(cmd, shell=True)
        
    # run hmmsearch
    def hmmsearch(self, for_out_reads, rev_out_reads, sequence_file_list, input_file_format, seq_type, threads, eval):
        suffix = [for_out_reads, rev_out_reads]
        table_title_list = []
    
        for seq_file in sequence_file_list:
            hmmout_table_title = suffix[0]
            table_title_list.append(hmmout_table_title)
            hmmsearch_cmd = " hmmsearch --cpu %s %s -o /dev/null --domtblout %s %s " % (threads, eval, hmmout_table_title, self.hmm)
            # TODO: capture stderr and report if the check_call fails
    
            if input_file_format == FORMAT_FASTA or input_file_format == FORMAT_FASTQ_GZ:
    
                if seq_type == 'P':
                    cmd = 'orfm %s | %s /dev/stdin' % (seq_file, hmmsearch_cmd)
                    # log
                    subprocess.check_call(["/bin/bash", "-c", cmd])
                elif seq_type == 'D':
    
                    if input_file_format == FORMAT_FASTQ_GZ:
                        cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s)) 2>&1 > /dev/null " % (hmmsearch_cmd, seq_file)
                        # log
                        subprocess.check_call(["/bin/bash", "-c", cmd])
    
                    elif input_file_format == FORMAT_FASTA:
                        cmd = "%s %s" % (hmmsearch_cmd, seq_file)
                        # log
                        subprocess.check_call(["/bin/bash", "-c", cmd])
    
                else:
                    Messenger().message('ERROR: Programming error')
                    exit(1)
    
    
            else:
                Messenger().message('ERROR: Suffix on %s not recegnised\n' % (seq_file))
                exit(1)
            del suffix[0]
    
        return table_title_list
        
    
    # run nhmmer
    def nhmmer(self, for_out_path, rev_out_path, sequence_file_list, input_file_format, threads, eval):    
        suffix = [for_out_path, rev_out_path]
        table_title_list = []
        for seq_file in sequence_file_list:
            hmmout_table_title = suffix.pop(0)
            table_title_list.append(hmmout_table_title)
            nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s %s" % (threads, eval, hmmout_table_title, self.hmm)
            
            if input_file_format == FORMAT_FASTA:
                cmd = "%s %s 2>&1 > /dev/null" % (nhmmer_cmd, seq_file)
                # log
    
                subprocess.check_call(["/bin/bash", "-c", cmd])
    
            elif input_file_format == FORMAT_FASTQ_GZ:
                cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s)) 2>&1 > /dev/null" % (nhmmer_cmd, seq_file)
                subprocess.check_call(["/bin/bash", "-c", cmd])
                # log
            else:
                Messenger().message('ERROR: Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (seq_file))
                exit(1)
    
        return table_title_list
    
    # Find Euk contmaination
    

    def check_euk_contamination(self, euk_free_path, out_table, reads, evals, avg_read_length, input_file_format, threads, eval, check_total_euks, raw_reads):
    
        
        contamination_list = []
        euk_uniq = []
        cutoff = 0.9*avg_read_length
        # do a nhmmer using a Euk specific hmm
        
        
        if check_total_euks:
            nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s /srv/db/graftm/0/HMM/Euk.hmm " % (threads, eval, out_table)
            
            if input_file_format == FORMAT_FASTA:
                cmd = nhmmer_cmd + raw_reads +' 2>&1 > /dev/null'
    
                # log
                subprocess.check_call(cmd, shell = True)
            
            elif input_file_format == FORMAT_FASTQ_GZ:
                cmd = nhmmer_cmd + " <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat " + raw_reads + ")) 2>&1 > /dev/null"
                # log
                subprocess.check_call(["/bin/bash", "-c", cmd])
        
    
        
            else:
                Messenger().error_message('Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (raw_reads))
                exit(1)
        
        else:    
            cmd = "nhmmer --cpu %s %s --tblout %s /srv/db/graftm/0/HMM/Euk.hmm %s 2>&1 > /dev/null " % (threads, eval, out_table, reads)
            # log
            subprocess.check_call(cmd, shell = True)
        
        
        
        # check for evalues that are lower, after eliminating hits with an alignment length of < 90% the length of the whole read.
        for line in open(out_table, 'r'):
            
            if not line.startswith('#'):
                
                line = line.split()
                try:
                    
                    if float(evals[line[0]][0]) >= float(line[12]):
                        
                        ali_length = float(line[6]) - float(line[7])
                        
                        if ali_length < 0:
                            ali_length = ali_length * -1.0
                        
                        if ali_length >= float(cutoff):
                            contamination_list.append(line[0])
                            euk_uniq.append(line[0])
                
                except KeyError:
    
                    if check_total_euks:
                        euk_uniq.append(line[0])
                    
                    else:
                        continue
                    
                    
        # Return Euk contamination
        if len(contamination_list) == 0:
            Messenger().message("No contaminating eukaryotic reads detected")
        
        else:
            Messenger().message("Found %s read(s) that may be eukaryotic, continuing without it/them" % len(contamination_list))
        
        # Write a file with the Euk free reads.
        with open(euk_free_path, 'w') as euk_free_output:
        
            for record in list(SeqIO.parse(open(reads, 'r'), 'fasta')):
    
                if record.id not in contamination_list:
                    SeqIO.write(record, euk_free_output, "fasta")
        if check_total_euks:
            euk_uniq = len(euk_uniq)
            contamination_list = len(contamination_list)
        else:
            euk_uniq = 'N/A'
            contamination_list = len(contamination_list)
            
        return contamination_list, euk_uniq

