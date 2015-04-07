import subprocess
import os
import json
import timeit

from Bio import SeqIO

from graftm.Messenger import Messenger
from graftm.assembler import TaxoGroup
from graftm.HouseKeeping import HouseKeeping


class Pplacer:
    ### Contains function related to processing alignment files to jplace files
    ### and running comparisons between forward and revere reads if reverse
    ### reads are provided.
    
    def __init__(self, refpkg):
        self.refpkg = refpkg
        self.HK = HouseKeeping()
    
    # Run pplacer
    def pplacer(self, output_file, output_path, input_path, threads, cmd_log):
        ## Runs pplacer on concatenated alignment file
        cmd = "pplacer -j %s --verbosity 0 --out-dir %s -c %s %s" % (threads, output_path, self.refpkg, input_path) # Set command
        self.HK.add_cmd(cmd_log, cmd) # Log it
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

    def guppy_class(self, main_guppy_path, jplace_list, cmd_log):
        ## Run guppy classify, and parse the output to the appropriate paths
        
        # Create concatenated guppy classify file from all .jplace files
        # created in the placement step
        cmd = 'guppy classify -c %s %s > %s' % (self.refpkg, ' '.join(jplace_list), main_guppy_path)
        self.HK.add_cmd(cmd_log, cmd)
        subprocess.check_call(cmd, shell=True)
        # Create list of guppys
        all_guppys = [x.rstrip() for x in open(main_guppy_path, 'r').readlines()]
        gup = []
        guppys = []
        for line in all_guppys:
            if 'name' in line and len(gup) == 0:
                gup.append(line)
            elif 'name' in line and len(gup) >= 0:
                guppys.append(gup)
                gup = [line]
            else:
                gup.append(line)
        guppys.append(gup)
        
        # Parse the guppy files.
        for idx, gup in enumerate(guppys):
            gup = [x for x in gup if x] # For each of the guppys remove empty components of the list
            out = os.path.join(os.path.dirname(jplace_list[idx]), 'placements.guppy') # Find the output
            r_num = len(list(set( [x.split()[0] for x in gup if 'name' not in x]))) # Calculate the number of placements
            with open(out, 'w') as out_guppy:
                for line in gup:
                    out_guppy.write(line + '\n')
        self.HK.delete([main_guppy_path])

        return
    
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
        
        stop = timeit.default_timer() # stop placement timer and log
        summary_dict['place_t'] = str( int(round((stop - start), 0)) )
        
        # Split the jplace file
        summary_dict = self.jplace_split(jplace, alias_hash, summary_dict)
        self.HK.delete([jplace])
        # Run guppy classify and parse the output 
        self.guppy_class(files.main_guppy_path(), summary_dict['jplace_path_list'], files.command_log_path())
        
        # If the reverse pipe has been specified, run the comparisons between the two pipelines. If not then just return.
        
        for base in summary_dict['base_list']:
            if summary_dict['reverse_pipe']:
                summary_dict[base] = Compare().compare_hits(summary_dict[base], base)
                summary_dict[base] = Compare().compare_placements(os.path.join(args.output_directory,base,'forward','placements.guppy'),
                                                                  os.path.join(args.output_directory,base,'reverse','placements.guppy'),
                                                                  summary_dict[base],
                                                                  args.placements_cutoff)
                return summary_dict
            
            elif not summary_dict['reverse_pipe']: # Set the trusted placements as 
                summary_dict[base]['trusted_placements'] = {}
                tc = TaxoGroup().guppy_splitter(os.path.join(args.output_directory,base,'placements.guppy'), args.placements_cutoff)
                for read, entry in tc.iteritems():
                    summary_dict[base]['trusted_placements'][read] = entry['placement']
        return summary_dict
        
class Compare:
    ### Functions for comparing forward and reverse read hits and placements
    
    def __init__(self): pass
    
    def compare_hits(self, hash, file_name):
        ## Take a paired read run, and compare hits between the two, report the 
        ## number of hits each, the crossover, and a 
        
        # Read in the read names for each file
        forward_read_names = set(hash['forward']['reads'].keys())
        reverse_read_names = set(hash['reverse']['reads'].keys())

        # Report and record the crossover
        crossover_hits = [x for x in forward_read_names if x in reverse_read_names]
        hash['crossover'] = crossover_hits
        Messenger().message("%s reads found that crossover in %s" % (str(len(crossover_hits)), 
                                                                   file_name))
        
        # Return the hash
        return hash
    
    def compare_placements(self, forward_guppy, reverse_guppy, hash, placement_cutoff):
        ## Take guppy placement file for the forward and reverse reads, compare
        ## the placement, and make a call as to which is to be trusted. Return
        ## a list of trusted reads for use by the summary step in GraftM
        
        # Read in guppy files
        forward_gup = TaxoGroup().guppy_splitter(forward_guppy, placement_cutoff)
        reverse_gup = TaxoGroup().guppy_splitter(reverse_guppy, placement_cutoff)
        
        # Report and record trusted placements
        comparison_hash = {'trusted_placements': {}} # Set up a hash that will record info on each placement
        for read in hash['crossover']:
            comparison_hash[read] = {} # make an entry for each read
            comparison_hash['trusted_placements'][read] = [] # Set up a taxonomy entry in trusted placements
            if len(forward_gup[read]['placement']) == len(reverse_gup[read]['placement']): # If the level of placement matches
                comparison_hash[read]['rank_length_match'] = True # Store True
            elif len(forward_gup[read]['placement']) != len(reverse_gup[read]['placement']):
                comparison_hash[read]['rank_length_match'] = False # Otherwise store False
            else:
                raise Exception('Programming Error: Comparison of placement resolution') 
            for idx, (f_rank, r_rank) in enumerate(zip(forward_gup[read]['placement'], reverse_gup[read]['placement'])): # For the each rank in the read placement
                if f_rank == r_rank: # If the classification at this rank matches
                    comparison_hash[read]['all_ranks_match'] = True # Store True
                    if comparison_hash[read]['rank_length_match']: # If both reads are to the same resolution
                        comparison_hash['trusted_placements'][read].append(f_rank) # Just append that rank
                    elif not comparison_hash[read]['rank_length_match']: # But if they are not at the same resolution
                        # check if we've reached the end of the forward or reverse taxonomy list
                        # and if we have, append the tail end of one with the longer taxonomy, because
                        # we can trust these placements are past the overall threshold, 
                        # but one has higher resolution and so is more useful to us.
                        if len(forward_gup[read]['placement'][idx:]) == 1: 
                            comparison_hash['trusted_placements'][read] += reverse_gup[read]['placement'][idx:]
                        elif len(reverse_gup[read]['placement'][idx:]) == 1:
                            comparison_hash['trusted_placements'][read] += forward_gup[read]['placement'][idx:]
                        elif len(forward_gup[read]['placement'][idx:]) > 1 and len(reverse_gup[read]['placement'][idx:]) > 1:
                            comparison_hash['trusted_placements'][read].append(f_rank)
                        else:
                            raise Exception('Programming Error')
                elif f_rank != r_rank: # Otherwise if the classification doesn't match
                    comparison_hash[read]['all_ranks_match'] = False #  store false
                    forward_confidence = forward_gup[read]['confidence'][idx] # Store confidence values for both directions
                    reverse_confidence = reverse_gup[read]['confidence'][idx]
                    if float(forward_confidence) > float(reverse_confidence): # If the forward read has more confidence
                        comparison_hash['trusted_placements'][read] += forward_gup[read]['placement'][idx:] # Store it
                    elif float(reverse_confidence) > float(forward_confidence): # Do the opposite if reverse read has more confidence
                        comparison_hash['trusted_placements'][read] += reverse_gup[read]['placement'][idx:] 
                    elif float(reverse_confidence) == float(forward_confidence): # If they are of the same value
                        break # Do nothing, because a decision cannot be made if the confidence is equal.
                    else:
                        raise Exception('Programming Error: Comparing confidence values')
                else:
                    raise Exception('Programming Error: Comparison of placement resolution')
        hash['comparison_hash'] = comparison_hash
        return hash # Return the hash
    
