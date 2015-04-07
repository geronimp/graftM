import subprocess
import os
import re
import timeit

from Bio import SeqIO
from collections import OrderedDict

from graftm.messenger import Messenger
from graftm.housekeeping import HouseKeeping

FORMAT_FASTA = 'FORMAT_FASTA'
FORMAT_FASTQ_GZ = 'FORMAT_FASTQ_GZ'

class Hmmer:

    def __init__(self, hmm):
        self.hmm = hmm
        self.hk = HouseKeeping()

    def hmmalign(self, input_path, run_stats, cmd_log, for_file, rev_file, for_sto_file, rev_sto_file, for_conv_file, rev_conv_file):
        # Align input reads to a specified hmm.
        if run_stats['rev_true']:
            read_info = run_stats['reads']
            reverse = []
            forward = []
            records = list(SeqIO.parse(open(input_path), 'fasta'))
            # Split the reads into reverse and forward lists
            for record in records:

                if read_info[record.id]['direction'] == '+':
                    forward.append(record)

                elif read_info[record.id]['direction'] == '-':
                    reverse.append(record)

                else:
                    raise Exception(Messenger().error_message('Programming error: hmmalign'))
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
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)
            cmd = 'hmmalign --trim -o %s %s %s 2>/dev/null; seqmagick convert %s %s' % (rev_sto_file,
                                                                                        self.hmm,
                                                                                        rev_file,
                                                                                        rev_sto_file,
                                                                                        rev_conv_file)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)

        # If there are only forward reads, just hmmalign and be done with it.
        else:
            cmd = 'hmmalign --trim -o %s %s %s ; seqmagick convert %s %s' % (for_sto_file,
                                                                             self.hmm,
                                                                             input_path,
                                                                             for_sto_file,
                                                                             for_conv_file)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)


    def hmmsearch(self, output_path, input_path, input_file_format, seq_type, threads, eval, cmd_log):
        ## Run a hmmsearch on the input_path raw reads, and return the name
        ## of the output table. Keep a log of the commands.
        # Define the base hmmsearch command.
        hmmsearch_cmd = "hmmsearch --cpu %s %s -o /dev/null --domtblout %s %s " % (threads, eval, output_path, self.hmm)
        # Choose an input to this base command based off the file format found.

        if seq_type == 'nucleotide': # If the input is nucleotide sequence
            cmd = 'orfm %s | %s /dev/stdin' % (input_path, hmmsearch_cmd) # call orfs on it, and search it
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])
        elif seq_type == 'protein': # If the input is amino acid sequence
            if input_file_format == FORMAT_FASTQ_GZ: # If its gzipped
                cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s)) 2>&1 > /dev/null " % (hmmsearch_cmd, input_path) # Unzip it and feed it into the base command
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])
            elif input_file_format == FORMAT_FASTA: # If it is in fasta format
                cmd = "%s %s" % (hmmsearch_cmd, input_path) # It can be searched directly, no manpulation required
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])
            else:
                raise Exception('Programming Error: error guessing input file format')
        else:
            raise Exception('Programming Error: error guessing input sequence type')
        return output_path

    def nhmmer(self, output_path, input_path, input_file_format, threads, eval, cmd_log):
        ## Run a nhmmer search on input_path file and return the name of
        ## resultant output table. Keep log of command.
        nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s %s" % (threads, eval, output_path, self.hmm)

        if input_file_format == FORMAT_FASTA:
            cmd = "%s %s 2>&1 > /dev/null" % (nhmmer_cmd, input_path)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])

        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s)) 2>&1 > /dev/null" % (nhmmer_cmd, input_path)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(["/bin/bash", "-c", cmd])

        else:
            raise Exception(Messenger().message('ERROR: Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (input_path)))

        return output_path

    # Filter Search
    def filter_nhmmer(self, file_path, threads, hmm_hash, cmd_log):
        suffix = [for_out_path, rev_out_path]
        table_title_list = []
        for seq_file in sequence_file_list:
            hmmout_table_title = suffix.pop(0)
            table_title_list.append(hmmout_table_title)
            nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s %s" % (threads, eval, hmmout_table_title, self.hmm)

            if input_file_format == FORMAT_FASTA:
                cmd = "%s %s 2>&1 > /dev/null" % (nhmmer_cmd, seq_file)
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])

            elif input_file_format == FORMAT_FASTQ_GZ:
                cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s)) 2>&1 > /dev/null" % (nhmmer_cmd, seq_file)
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])

            else:
                Messenger().message('ERROR: Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (seq_file))
                exit(1)

        return table_title_list
    # Find Euk contmaination


    def check_euk_contamination(self, output_path, euk_free_output_path, input_path, run_stats, input_file_format, check_total_euks, threads, evalue, raw_reads, base, cmd_log, euk_hmm):
        reads_with_better_euk_hit = []
        reads_unique_to_eukaryotes = []
        cutoff = 0.9*run_stats['read_length']
        # do a nhmmer using a Euk specific hmm
        if check_total_euks:
            nhmmer_cmd = "nhmmer --cpu %s %s --tblout %s %s" % (threads, evalue, output_path, euk_hmm)

            if input_file_format == FORMAT_FASTA:
                cmd = "%s %s 2>&1 > /dev/null" % (nhmmer_cmd, raw_reads)
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(cmd, shell = True)

            elif input_file_format == FORMAT_FASTQ_GZ:
                cmd = "%s <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat %s )) 2>&1 > /dev/null" %  (nhmmer_cmd, raw_reads)
                self.hk.add_cmd(cmd_log, cmd)
                subprocess.check_call(["/bin/bash", "-c", cmd])

            else:
                raise Exception(Messenger().error_message('Suffix on %s not familiar. Please submit an .fq.gz or .fa file\n' % (raw_reads)))

        else:
            cmd = "nhmmer --cpu %s %s --tblout %s /srv/db/graftm/0/HMM/Euk.hmm %s 2>&1 > /dev/null " % (threads, evalue, output_path, input_path)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell = True)

        # check for evalues that are lower, after eliminating hits with an
        # alignment length of < 90% the length of the whole read.
        for line in [x.rstrip().split() for x in open(output_path, 'r') if not x.startswith('#')]:
            read_name = line[0]
            eval = line[12]
            ali_length = float(line[6]) - float(line[7])

            try:
                if float(run_stats['reads'][read_name]['evalue']) >= float(eval):
                    if ali_length < 0:
                        ali_length = ali_length * -1.0

                    if ali_length >= float(cutoff):
                        reads_with_better_euk_hit.append(read_name)
                        reads_unique_to_eukaryotes.append(read_name)

            except KeyError:
                if check_total_euks:
                    reads_unique_to_eukaryotes.append(read_name)

                else:
                    continue

        # Return Euk contamination
        if len(reads_unique_to_eukaryotes) == 0:
            Messenger().message("No contaminating eukaryotic reads detected in %s" % (os.path.basename(raw_reads)))

        else:
            Messenger().message("Found %s read(s) that may be eukaryotic, continuing without it/them" % len(reads_unique_to_eukaryotes))

        # Write a file with the Euk free reads.
        with open(euk_free_output_path, 'w') as euk_free_output:
            for record in list(SeqIO.parse(open(input_path, 'r'), 'fasta')):
                if record.id not in reads_unique_to_eukaryotes:
                    SeqIO.write(record, euk_free_output, "fasta")

        if check_total_euks:
            run_stats['euk_uniq'] = len(euk_uniq)
            run_stats['euk_contamination'] = len(reads_unique_to_eukaryotes)

        else:
            run_stats['euk_uniq'] = 'N/A'
            run_stats['euk_contamination'] = len(reads_unique_to_eukaryotes)

        return run_stats, euk_free_output_path

    def filter_hmmsearch(self, output_hash, contents, args, input_file_format, cmd_log):
        for seq_file in sequence_file_list:
            hmmout_table_title = suffix[0]
            table_title_list.append(hmmout_table_title)
            hmmsearch_cmd = " hmmsearch --cpu %s %s -o /dev/null --domtblout %s %s " % (threads, eval, hmmout_table_title, self.hmm)
            # TODO: capture stderr and report if the check_call fails

            if input_file_format == FORMAT_FASTA or input_file_format == FORMAT_FASTQ_GZ:

                if contents.pipe == 'P':
                    cmd = 'orfm %s | %s /dev/stdin' % (seq_file, hmmsearch_cmd)
                    self.hk.add_cmd(cmd_log, cmd)
                    subprocess.check_call(["/bin/bash", "-c", cmd])



            else:
                Messenger().message('ERROR: Suffix on %s not recegnised\n' % (seq_file))
                exit(1)
            del suffix[0]

        return table_title_list

    def csv_to_titles(self, output_path, input_path, graftm_pipeline, run_stats, seq_type):
        ## process hmmsearch/nhmmer results into a list of titles to <base_filename>_readnames.txt
        # Assume there are no reverse reads detected
        run_stats['rev_true'] = False
        reads_list = []
        orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
        run_stats['reads'] = {}
        hmm_table = [x.rstrip().split() for x in open(input_path, 'r').readlines() if not x.startswith('#')]
        for entry in hmm_table:
            read_name = entry[0]
            evalue = entry[12]
            direction = entry[11]
            if direction == '-':
                run_stats['rev_true'] = True
            if graftm_pipeline == 'D':
                run_stats['reads'][read_name] = {'evalue': evalue,
                                                 'direction': direction}
                reads_list.append(read_name)
            elif graftm_pipeline == 'P':
                if seq_type == 'nucleotide':
                # The original reads file contains sequences like
                # >eg and comment
                # where orfm gives orfs in the form of
                # >eg_1_2_3 and comment
                # The read_name here is the orfm style, we want to add to the titles_list
                # the original form
                    regex_match = orfm_regex.match(read_name)
                    if regex_match is not None:
                        read_name = regex_match.groups(0)[0]
                        run_stats['reads'][read_name] = {'evalue': evalue,
                                                         'direction': direction}
                        reads_list.append(read_name)
                    else:
                        raise Exception("Unexpected form of ORF name found: %s" % read_name)
                else:
                    reads_list.append(read_name)

        if len(reads_list) > 0:
            Messenger().message('Found %s read(s) in %s' % (len(reads_list), os.path.basename(input_path)))
        else:
            Messenger().message('%s reads found, cannot continue' % (len(reads_list)))
            raise Exception('No reads to continue with')

        ## write the read names to output
        with open(output_path, 'w') as output_file:
            for record in reads_list:
                output_file.write(record + '\n')

        return run_stats, output_path

    def extract_from_raw_reads(self, output_path, input_path, raw_sequences_path, input_file_format, cmd_log):
        # Use the readnames specified to extract from the original sequence
        # file to a fasta formatted file.

        # Run fxtract to obtain reads form original sequence file
        fxtract_cmd = "fxtract -H -X -f %s " % input_path

        if input_file_format == FORMAT_FASTA:
            cmd = "%s %s > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)

        elif input_file_format == FORMAT_FASTQ_GZ:
            cmd = "%s %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > %s" % (fxtract_cmd, raw_sequences_path, output_path)
            self.hk.add_cmd(cmd_log, cmd)
            subprocess.check_call(cmd, shell=True)

        else:
            raise Exception("Programming error")

        return output_path

    def check_read_length(self, reads, pipe):
        lengths = []
        record_list = []
        suffixes = ['_for', '_rev']
        # First check if the reverse pipe is happening, because the read names
        # are different.

        record_list += list(SeqIO.parse(open(reads, 'r'), 'fasta'))

        for record in record_list:
            lengths.append(len(record.seq))
        if pipe == "P":
            return (sum(lengths) / float(len(lengths)))/3
        elif pipe =="D":
            return sum(lengths) / float(len(lengths))


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
                orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
                regex_match = orfm_regex.match(fasta_id)

                if regex_match is not None:
                    output_file.write(regex_match.groups(0)[0] + '\n')
                else:
                    output_file.write(fasta_id)

                output_file.write(fasta_seq)

    def extract_orfs(self, input_path, raw_orf_path, hmmsearch_out_path, orf_titles_path, orf_out_path, hmm, cmd_log):
        ## Extract only the orfs that hit the hmm, return sequence file with
        ## within.

        # Call orfs on the sequences
        cmd = 'orfm %s > %s' % (input_path, raw_orf_path)
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)

        # Search for the correct reading fram
        cmd = 'hmmsearch -o /dev/null --tblout %s %s %s' % (hmmsearch_out_path, hmm, raw_orf_path)
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)

        # Extract hit orf names
        with open(orf_titles_path, 'w') as output:
            for title in [x.split(' ')[0] for x in open(hmmsearch_out_path).readlines() if not x.startswith('#')]:
                output.write(str(title) + '\n')

        # Extract the reads using the titles.
        cmd = 'fxtract -H -X -f %s %s > %s' % (orf_titles_path, raw_orf_path, orf_out_path)
        self.hk.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)

        # Return name of output file
        return orf_out_path

    def p_search(self, files, args, run_stats, base, input_file_format, raw_reads):
        # Main pipe of search step in protein pipeline:
        # recieves reads, and returns hits
        start = timeit.default_timer() # Start search timer
        # Searching raw reads with HMM
        hit_table = self.hmmsearch(files.hmmsearch_output_path(base),
                                   raw_reads,
                                   input_file_format,
                                   args.input_sequence_type,
                                   args.threads,
                                   args.eval,
                                   files.command_log_path())
        # Processing the output table to give you the readnames of the hits
        run_stats, hit_readnames= self.csv_to_titles(files.readnames_output_path(base),
                                                     hit_table,
                                                     args.type,
                                                     run_stats,
                                                     args.input_sequence_type)
        # Extract the hits form the original raw read file
        hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                hit_readnames,
                                                raw_reads,
                                                input_file_format,
                                                files.command_log_path())
        # Extract the orfs of these reads that hit the original search
        hit_orfs = self.extract_orfs(hit_reads,
                                     files.orf_output_path(base),
                                     files.orf_hmmsearch_output_path(base),
                                     files.orf_titles_output_path(base),
                                     files.orf_fasta_output_path(base),
                                     args.hmm_file,
                                     files.command_log_path())
        # Define the average read length of the hits
        run_stats['read_length'] = self.check_read_length(hit_orfs, "P")
        # Stop and log search timer
        stop = timeit.default_timer()
        run_stats['search_t'] = str(int(round((stop - start), 0)) )
        # Return hit reads, and summary hash
        return hit_orfs, run_stats

    def d_search(self, files, args, run_stats, base, input_file_format, raw_reads):
        # Main pipe of search step in nucleotide pipeline:
        # recieves reads, and returns hits
        start = timeit.default_timer() # Start search timer

        # First search the reads using the HMM
        hit_table = self.nhmmer(files.hmmsearch_output_path(base),
                                raw_reads,
                                input_file_format,
                                args.threads,
                                args.eval,
                                files.command_log_path())

        # Next, get a list of readnames
        run_stats, hit_readnames= self.csv_to_titles(files.readnames_output_path(base),
                                                     hit_table,
                                                     args.type,
                                                     run_stats,
                                                     None)

        # And extract them from the original sequence file
        hit_reads = self.extract_from_raw_reads(files.fa_output_path(base),
                                                hit_readnames,
                                                raw_reads,
                                                input_file_format,
                                                files.command_log_path())

        # Define the read length
        run_stats['read_length'] = self.check_read_length(hit_reads, "D")

        # Stop timing search and start timing euk check step.
        stop = timeit.default_timer()
        run_stats['search_t'] = str(int(round((stop - start), 0)) )
        start = timeit.default_timer()

        # Check for Eukarytoic contamination
        Messenger().message("Checking for Eukaryotic contamination")
        run_stats, hit_reads = self.check_euk_contamination(files.euk_contam_path(base),
                                                            files.euk_free_path(base),
                                                            hit_reads,
                                                            run_stats,
                                                            input_file_format,
                                                            args.check_total_euks,
                                                            args.threads,
                                                            args.eval,
                                                            raw_reads,
                                                            base,
                                                            files.command_log_path(),
                                                            args.euk_hmm_file)

        # Stop timing eukaryote check
        stop = timeit.default_timer()
        run_stats['euk_check_t'] = str(int(round((stop - start), 0)) )

        # Finally, return the hits
        return hit_reads, run_stats

    def align(self, files, args, run_stats, base, reads):
        # This pipeline takes unaligned reads, and aligns them agains a hmm,
        # regardless of their direction. Aligned reads with base insertions
        # removed are returned in the end. Times and commands are logged.

        start = timeit.default_timer()

        # HMMalign the forward reads, and reverse complement reads.
        self.hmmalign(reads,
                      run_stats,
                      files.command_log_path(),
                      files.output_for_path(base),
                      files.output_rev_path(base),
                      files.sto_for_output_path(base),
                      files.sto_rev_output_path(base),
                      files.conv_output_for_path(base),
                      files.conv_output_rev_path(base))

        # Correct the alignment for base insertions.
        if run_stats['rev_true']:
            self.alignment_correcter([files.conv_output_for_path(base), files.conv_output_rev_path(base)],
                                     files.aligned_fasta_output_path(base))
        else:
            self.alignment_correcter([files.conv_output_for_path(base)],
                                     files.aligned_fasta_output_path(base))
        stop = timeit.default_timer()
        run_stats['aln_t'] = str(int(round((stop - start), 0)) )

        # Return
        return files.aligned_fasta_output_path(base), run_stats

