import subprocess
import os
import json
import timeit

from Bio import SeqIO

from graftm.classify import Classify
from graftm.messenger import Messenger
from graftm.assembler import TaxoGroup
from graftm.housekeeping import HouseKeeping


class Pplacer:
    ### Contains function related to processing alignment files to jplace files
    ### and running comparisons between forward and revere reads if reverse
    ### reads are provided.

    def __init__(self, refpkg):
        self.refpkg = refpkg
        self.hk = HouseKeeping()

    # Run pplacer
    def pplacer(self, output_file, output_path, input_path, threads, cmd_log):
        ## Runs pplacer on concatenated alignment file
        cmd = "pplacer -j %s --verbosity 0 --out-dir %s -c %s %s" % (threads, output_path, self.refpkg, input_path) # Set command
        self.hk.add_cmd(cmd_log, cmd) # Log it
        subprocess.check_call(cmd, shell=True) # Run it
        output_path = '.'.join(input_path.split('.')[:-1]) + '.jplace'
        return output_path

    def alignment_merger(self, alignment_files, output_alignment_path):
        ## Concatenate aligned read_files into one file. Each read with it's
        ## own unique identifier assigning it to a particular origin file

        alias_hash = {} # Set up a hash with file names and their unique identifier
        file_number = 0 # file counter (unique identifier)
        with open(output_alignment_path, 'w') as output:
            for alignment_file in alignment_files: # For each alignment
                alignments = list(SeqIO.parse(open(alignment_file, 'r'), 'fasta')) # read list
                for record in alignments: # For each record in the read list
                    record.id = record.id + '_' + str(file_number) # append the unique identifier to the record id
                SeqIO.write(alignments, output, "fasta") # And write the reads to the file
                alias_hash[str(file_number)] = {'output_path': os.path.join(os.path.dirname(alignment_file),'placements.jplace') ,
                                             'place': []}
                file_number += 1
        return alias_hash

    def jplace_split(self, jplace_file, alias_hash, summary_dict):
        ## Split the jplace file into their respective directories

        # Load the placement file
        placement_file = json.load(open(jplace_file))

        # Parse the placements based on unique identifies appended to the end
        # of each read
        for placement in placement_file['placements']: # for each placement
            hash = {} # create an empty hash
            for alias in alias_hash: # For each alias, append to the 'place' list each read that identifier
                hash = {'p': placement['p'],
                        'nm': [nm for nm in placement['nm'] if nm[0].split('_')[-1] == alias]}
                alias_hash[alias]['place'].append(hash)

        # Write the jplace file to their respective file paths.
        jplace_path_list = []
        for alias in alias_hash:
            output = {'fields': placement_file['fields'],
                      'version': placement_file['version'],
                      'tree':  placement_file['tree'],
                      'placements': alias_hash[alias]['place'],
                      'metadata': placement_file['metadata']}
            with open(alias_hash[alias]['output_path'], 'w') as output_path:
                json.dump(output, output_path, ensure_ascii=False)
            jplace_path_list.append(alias_hash[alias]['output_path'])
        summary_dict['jplace_path_list'] = jplace_path_list
        return summary_dict

    def place(self, summary_dict, files, args):
        ## Pipeline taking multiple alignment files and returning multiple
        ## placement and guppy files, as well as the comparison between forward
        ## and reverse reads, if the reverse pipeline is selected

        start = timeit.default_timer() # Start placement timer

        # Merge the alignments so they can all be placed at once.
        alias_hash = self.alignment_merger(summary_dict['seqs_list'], files.comb_aln_fa())
        
        # Run pplacer on merged file
        jplace = self.pplacer(files.jplace_output_path(), args.output_directory, files.comb_aln_fa(), args.threads, files.command_log_path())
        Messenger().message("Placements finished")

        stop = timeit.default_timer() # stop placement timer and log
        summary_dict['place_t'] = str( int(round((stop - start), 0)) )

        # Split the jplace file
        summary_dict = self.jplace_split(jplace, alias_hash, summary_dict)
        
        #Read the json of refpkg
        Messenger().message("Reading classifications")
        tax_descr=json.load(open(self.refpkg+'/CONTENTS.json'))['files']['taxonomy']
        classifications=Classify(os.path.join(self.refpkg,tax_descr)).assignPlacement(jplace, args.placements_cutoff, 'reads')
        Messenger().message("Reads classified.")
        
        # If the reverse pipe has been specified, run the comparisons between the two pipelines. If not then just return.

        for idx, base in enumerate(summary_dict['base_list']):
            if summary_dict['reverse_pipe']:
                summary_dict[base] = Compare().compare_hits(summary_dict[base], base)

                forward_gup=classifications.pop(sorted(classifications.keys())[0]) 
                reverse_gup=classifications.pop(sorted(classifications.keys())[0])
                summary_dict[base] = Compare().compare_placements(forward_gup,
                                                                  reverse_gup,
                                                                  summary_dict[base],
                                                                  args.placements_cutoff)

            elif not summary_dict['reverse_pipe']: # Set the trusted placements as
                summary_dict[base]['trusted_placements'] = {}
                for read, entry in classifications[str(idx)].iteritems():
                    summary_dict[base]['trusted_placements'][read] = entry['placement']
        return summary_dict

class Compare:
    ### Functions for comparing forward and reverse read hits and placements

    def __init__(self): pass

    def compare_hits(self, hash, file_name):
        ## Take a paired read run, and compare hits between the two, report the
        ## number of hits each, the crossover, and a
        def check_reads(read_lists):
            for read_list in read_lists:
                read_check = [read for read in read_list if read.startswith('FCC')]
                if read_check:
                    return True
                elif not read_check:
                    continue
            return False
        # Read in the read names for each file
        check = check_reads([hash['forward']['reads'].keys(), hash['reverse']['reads'].keys()])
        if check:
            forward_read_names = [x.replace('/1', '') for x in hash['forward']['reads'].keys()]
            reverse_read_names = [x.replace('/2', '') for x in hash['reverse']['reads'].keys()]
        elif not check:
            forward_read_names = set(hash['forward']['reads'].keys())
            reverse_read_names = set(hash['reverse']['reads'].keys())
        # Report and record the crossover
        crossover_hits = [x for x in forward_read_names if x in reverse_read_names]
        hash['crossover'] = crossover_hits
        # Check if there are reads to continue with.
        if len(crossover_hits) > 0:
            Messenger().message("%s reads found that crossover in %s" % (str(len(crossover_hits)),
                                                                       file_name))
        elif len(crossover_hits) == 0:
            Messenger().message("%s reads found that crossover in %s, No reads to use!" % (str(len(crossover_hits)),
                                                                                           file_name))            
        else:
            raise Exception
        # Return the hash

        return hash 

    def compare_placements(self, forward_gup, reverse_gup, hash, placement_cutoff):
        ## Take guppy placement file for the forward and reverse reads, compare
        ## the placement, and make a call as to which is to be trusted. Return
        ## a list of trusted reads for use by the summary step in GraftM
        # Report and record trusted placements
        comparison_hash = {'trusted_placements': {}} # Set up a hash that will record info on each placement
        for read in hash['crossover']:
            if read.startswith('FCC'):
                f_read = read + '/1'
                r_read = read + '/2'
            else:
                f_read = read
                r_read = read
            # Check read was placed
            if forward_gup.get(f_read) is None or reverse_gup.get(r_read) is None:
                Messenger().message('Warning: %s was not inserted into tree' % str(f_read))
                continue # Skip read for now
            comparison_hash[read] = {} # make an entry for each read
            comparison_hash['trusted_placements'][read] = [] # Set up a taxonomy entry in trusted placements

            if len(forward_gup[f_read]['placement']) == len(reverse_gup[r_read]['placement']): # If the level of placement matches
                comparison_hash[read]['rank_length_match'] = True # Store True
            elif len(forward_gup[f_read]['placement']) != len(reverse_gup[r_read]['placement']):
                comparison_hash[read]['rank_length_match'] = False # Otherwise store False
            else:
                raise Exception('Programming Error: Comparison of placement resolution')
            for idx, (f_rank, r_rank) in enumerate(zip(forward_gup[f_read]['placement'], reverse_gup[r_read]['placement'])): # For the each rank in the read placement                    
                if f_rank == r_rank: # If the classification at this rank matches
                    comparison_hash[read]['all_ranks_match'] = True # Maintain the all ranks match are true
                    if comparison_hash[read]['rank_length_match']: # If both reads are to the same resolution
                        comparison_hash['trusted_placements'][read].append(f_rank) # Just append that rank
                    elif not comparison_hash[read]['rank_length_match']: # But if they are not at the same resolution
                        # check if we've reached the end of the forward or reverse taxonomy list
                        # and if we have, append the tail end of one with the longer taxonomy, because
                        # we can trust these placements are past the overall threshold,
                        # but one has higher resolution and so is more useful to us.
                        if len(forward_gup[f_read]['placement'][idx:]) == 1:
                            comparison_hash['trusted_placements'][read] += reverse_gup[r_read]['placement'][idx:]
                        elif len(reverse_gup[r_read]['placement'][idx:]) == 1:
                            comparison_hash['trusted_placements'][read] += forward_gup[f_read]['placement'][idx:]
                        elif len(forward_gup[f_read]['placement'][idx:]) > 1 and len(reverse_gup[r_read]['placement'][idx:]) > 1:
                            comparison_hash['trusted_placements'][read].append(f_rank)
                        else:
                            raise Exception('Programming Error')
                elif f_rank != r_rank: # Otherwise if the classification doesn't match
                    comparison_hash[read]['all_ranks_match'] = False # Store that all ranks do not match up
                    forward_confidence = forward_gup[f_read]['confidence'][idx] # Store confidence values for both directions
                    reverse_confidence = reverse_gup[r_read]['confidence'][idx]
                    if float(forward_confidence) > float(reverse_confidence): # If the forward read has more confidence
                        comparison_hash['trusted_placements'][read] += forward_gup[f_read]['placement'][idx:] # Store the taxonomy of that read from that point on
                        break
                    elif float(reverse_confidence) > float(forward_confidence): # Do the opposite if reverse read has more confidence
                        comparison_hash['trusted_placements'][read] += reverse_gup[r_read]['placement'][idx:]
                        break
                    elif float(reverse_confidence) == float(forward_confidence): # If they are of the same value
                        break # Do nothing, because a decision cannot be made if the confidence is equal.
                    else:
                        raise Exception('Programming Error: Comparing confidence values')
                else:
                    raise Exception('Programming Error: Comparison of placement resolution')               

        hash['comparison_hash'] = comparison_hash
        return hash # Return the hash

