import subprocess
import os
import json
import logging
import time
import re

from Bio import SeqIO

from graftm.timeit import Timer
from graftm.classify import Classify
from graftm.housekeeping import HouseKeeping
T=Timer()


class Pplacer:
    ### Contains function related to processing alignment files to jplace files
    ### and running comparisons between forward and revere reads if reverse
    ### reads are provided.

    def __init__(self, refpkg):
        self.refpkg = refpkg
        self.hk = HouseKeeping()
        
    # Run pplacer
    def pplacer(self, output_file, output_path, input_path, threads):
        ## Runs pplacer on concatenated alignment file
        cmd = "pplacer -j %s --verbosity 0 --out-dir %s -c %s %s" % (str(threads), output_path, self.refpkg, input_path) # Set command
        logging.debug("Running command: %s" % cmd)
        # Log it
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
                alias_hash[str(file_number)] = {'output_path': os.path.join(os.path.dirname(alignment_file),'placements.jplace'),
                                                'place': []}
                file_number += 1
        return alias_hash

    def jplace_split(self, jplace_file, alias_hash):
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
            output = {'fields'     : placement_file['fields'],
                      'version'    : placement_file['version'],
                      'tree'       : placement_file['tree'],
                      'placements' : alias_hash[alias]['place'],
                      'metadata'   : placement_file['metadata']}
            with open(alias_hash[alias]['output_path'], 'w') as output_path:
                json.dump(output, output_path, ensure_ascii=False)
            jplace_path_list.append(alias_hash[alias]['output_path'])
        return jplace_path_list
    
    @T.timeit
    def place(self, reverse_pipe, seqs_list, resolve_placements, files, args,
              slash_endings, tax_descr):
        '''
        placement - This is the placement pipeline in GraftM, in aligned reads 
                    are placed into phylogenetic trees, and the results interpreted.
                    If reverse reads are used, this is where the comparisons are made
                    between placements, for the summary tables to be build in the
                    next stage.
         
        Parameters
        ----------
        reverse_pipe : bool
            True: reverse reads are placed separately
            False: no reverse reads to place.
        seqs_list : list
            list of paths to alignment fastas to be placed into the tree
        resolve_placements : bool
            True:resolve placements to their most trusted taxonomy
            False: classify reads to their most trusted taxonomy, until the 
                   confidence cutoff is reached. 
        files : list
            graftM output file name object
        args : obj
            argparse object
        Returns
        ------- 
        trusted_placements : dict
            dictionary of reads and their trusted placements
        '''
        trusted_placements = {}
        
        # Merge the alignments so they can all be placed at once.
        alias_hash = self.alignment_merger(seqs_list, files.comb_aln_fa())
        
        # Run pplacer on merged file
        jplace = self.pplacer(files.jplace_output_path(), args.output_directory, files.comb_aln_fa(), args.threads)
        logging.info("Placements finished")

        # Split the jplace file
        self.jplace_split(jplace, alias_hash)
        
        #Read the json of refpkg
        logging.info("Reading classifications")
        classifications=Classify(tax_descr).assignPlacement(
                                                           jplace, 
                                                           args.placements_cutoff, 
                                                           'reads', 
                                                           resolve_placements
                                                           )
                                         
        self.hk.delete([jplace])# Remove combined split, not really useful
        logging.info("Reads classified.")
        # If the reverse pipe has been specified, run the comparisons between the two pipelines. If not then just return.
        
        for idx, file in enumerate(seqs_list):
            
            if reverse_pipe:
                base_file=os.path.basename(file).replace('_forward_hits.aln.fa', '')
                forward_gup=classifications.pop(sorted(classifications.keys())[0]) 
                reverse_gup=classifications.pop(sorted(classifications.keys())[0])
                seqs_list.pop(idx+1)
                placements_hash = Compare().compare_placements(
                                                               forward_gup,
                                                               reverse_gup,
                                                               args.placements_cutoff,
                                                               slash_endings,
                                                               base_file
                                                               )
                trusted_placements[base_file]=placements_hash['trusted_placements']
                
            else: # Set the trusted placements as
                base_file=os.path.basename(file).replace('_hits.aln.fa', '')
                trusted_placements[base_file]={}
                for read, entry in classifications[str(idx)].iteritems():
                    trusted_placements[base_file][read] = entry['placement']
        
        return trusted_placements

class Compare:
    ### Functions for comparing forward and reverse read hits and placements

    def __init__(self): pass

    def _compare_hits(self, forward_reads, reverse_reads, file_name, slash_endings):
        ## Take a paired read run, and compare hits between the two, report the
        ## number of hits each, the crossover, and a
        
        def remove_endings(read_list, slash_endings):
            orfm_regex = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
            d = {}
            for read in read_list:
                orfm_match=orfm_regex.match(read)
                
                if orfm_match:
                    if slash_endings:
                        new_read=orfm_match.groups(0)[0][:-2]
                    else:
                        new_read=orfm_match.groups(0)[0]
                elif slash_endings:
                    new_read=read[:-2]
                else:
                    new_read=read
                d[new_read]=read
            return d
        
        forward_reads=remove_endings(forward_reads, slash_endings)
        reverse_reads=remove_endings(reverse_reads, slash_endings)
        # Report and record the crossover
        crossover_hits = [x for x in forward_reads.keys() if x in reverse_reads.keys()]

        # Check if there are reads to continue with.
        if len(crossover_hits) > 0:
            logging.info("%s reads found that crossover in %s" % (str(len(crossover_hits)),
                                                                       file_name))
        elif len(crossover_hits) == 0:
            logging.info("%s reads found that crossover in %s, No reads to use!" % (str(len(crossover_hits)),
                                                                                           file_name))            
        else:
            raise Exception
        # Return the hash

        return crossover_hits, forward_reads, reverse_reads

    def compare_placements(self, forward_gup, reverse_gup, placement_cutoff, slash_endings, base_file):
        ## Take guppy placement file for the forward and reverse reads, compare
        ## the placement, and make a call as to which is to be trusted. Return
        ## a list of trusted reads for use by the summary step in GraftM
        # Report and record trusted placements        
        crossover, for_dict, rev_dict = self._compare_hits(forward_gup.keys(), 
                                                           reverse_gup.keys(), 
                                                           base_file, 
                                                           slash_endings)
        comparison_hash = {'trusted_placements': {}} # Set up a hash that will record info on each placement
        for read in crossover:
            f_read = for_dict[read]
            r_read = rev_dict[read]

            # Check read was placed
            if forward_gup.get(f_read) is None or reverse_gup.get(r_read) is None:
                logging.info('Warning: %s was not inserted into tree' % str(f_read))
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
                
        return comparison_hash # Return the hash

