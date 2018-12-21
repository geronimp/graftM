import extern
import os
import json
import logging
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
        extern.run(cmd)
        output_path = '.'.join(input_path.split('.')[:-1]) + '.jplace'
        return output_path

    def alignment_merger(self, alignment_files, output_alignment_path):
        ## Concatenate aligned read_files into one file. Each read with it's
        ## own unique identifier assigning it to a particular origin file
        alias_hash = {} # Set up a hash with file names and their unique identifier
        file_number = 0 # file counter (unique identifier)
        with open(output_alignment_path, 'w') as output:
            for alignment_file in alignment_files: # For each alignment
                if alignment_file is not None:
                    alignments = list(SeqIO.parse(open(alignment_file, 'r'), 'fasta')) # read list
                    for record in alignments: # For each record in the read list
                        record.id = record.id + '_' + str(file_number) # append the unique identifier to the record id
                    SeqIO.write(alignments, output, "fasta") # And write the reads to the file
                    alias_hash[str(file_number)] = {'output_path': os.path.join(os.path.dirname(alignment_file), 'placements.jplace')}
                file_number += 1
        return alias_hash

    def convert_cluster_dict_keys_to_aliases(self, cluster_dict, alias_hash):
        '''
        Parameters
        ----------
        cluster_dict : dict
            dictionary stores information on pre-placement clustering
        alias_hash : dict
            Stores information on each input read file given to GraftM, the
            corresponding reads found within each file, and their taxonomy

        Returns
        --------
        updated cluster_dict dict containing alias indexes for keys.
        '''
        output_dict = {}
        directory_to_index_dict = {os.path.split(item["output_path"])[0] : key
                            for key, item in alias_hash.iteritems()}
        for key, item in cluster_dict.iteritems():
            cluster_file_directory = os.path.split(key)[0]
            cluster_idx = directory_to_index_dict[cluster_file_directory]
            output_dict[cluster_idx] = item
        return output_dict



    def jplace_split(self, original_jplace, cluster_dict):
        '''
        To make GraftM more efficient, reads are dereplicated and merged into
        one file prior to placement using pplacer. This function separates the
        single jplace file produced by this process into the separate jplace
        files, one per input file (if multiple were provided) and backfills
        abundance (re-replicates?) into the placement file so analyses can be
        done using the placement files.

        Parameters
        ----------
        original_jplace : dict (json)
            json .jplace file from the pplacer step.
        cluster_dict : dict
            dictionary stores information on pre-placement clustering

        Returns
        -------
        A dict containing placement hashes to write to
        new jplace file. Each key represents a file alias
        '''
        output_hash = {}

        for placement in original_jplace['placements']: # for each placement
            alias_placements_list = []
            nm_dict = {}

            p = placement['p']
            if 'nm' in placement.keys():
                nm = placement['nm']
            elif 'n' in placement.keys():
                nm = placement['n']
            else:
                raise Exception("Unexpected jplace format: Either 'nm' or 'n' are expected as keys in placement jplace .JSON file")

            for nm_entry in nm:
                nm_list = []
                placement_read_name, plval = nm_entry
                read_alias_idx = placement_read_name.split('_')[-1] # Split the alias
                                    # index out of the read name, which
                                    # corresponds to the input file from
                                    # which the read originated.
                read_name = '_'.join(placement_read_name.split('_')[:-1])
                read_cluster = cluster_dict[read_alias_idx][read_name]
                for read in read_cluster:
                    nm_list.append([read.name, plval])
                if read_alias_idx not in nm_dict:
                    nm_dict[read_alias_idx] = nm_list
                else:
                    nm_dict[read_alias_idx] += nm_entry

            for alias_idx, nm_list in nm_dict.iteritems():
                placement_hash = {'p': p,
                                  'nm': nm_list}
                if alias_idx not in output_hash:
                    output_hash[alias_idx] = [placement_hash]
                else:
                    output_hash[alias_idx].append(placement_hash)
        return output_hash

    def write_jplace(self, original_jplace, alias_hash):
        # Write the jplace file to their respective file paths.
        for alias_idx in alias_hash.keys():
            output = {'fields'     : original_jplace['fields'],
                      'version'    : original_jplace['version'],
                      'tree'       : original_jplace['tree'],
                      'placements' : alias_hash[alias_idx]['place'],
                      'metadata'   : original_jplace['metadata']}
            with open(alias_hash[alias_idx]['output_path'], 'w') as output_io:
                json.dump(output, output_io, ensure_ascii=False, indent=3, separators=(',', ': '))

    @T.timeit
    def place(self, reverse_pipe, seqs_list, resolve_placements, files, args,
              slash_endings, tax_descr, clusterer):
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
        files_to_delete = []
        # Merge the alignments so they can all be placed at once.
        alias_hash = self.alignment_merger(seqs_list, files.comb_aln_fa())
        files_to_delete += seqs_list
        files_to_delete.append(files.comb_aln_fa())
        if os.path.getsize(files.comb_aln_fa()) == 0:
            logging.debug("Combined alignment file has 0 size, not running pplacer")
            to_return = {}
            for idx, file in enumerate(seqs_list):
                base_file=os.path.basename(file).replace('_forward_hits.aln.fa', '')
                to_return[base_file] = {}
            return to_return

        # Run pplacer on merged file
        jplace = self.pplacer(files.jplace_output_path(), args.output_directory, files.comb_aln_fa(), args.threads)
        files_to_delete.append(jplace)
        logging.info("Placements finished")

        #Read the json of refpkg
        logging.info("Reading classifications")
        classifications=Classify(tax_descr).assignPlacement(
                                                           jplace,
                                                           args.placements_cutoff,
                                                           resolve_placements
                                                           )
        logging.info("Reads classified")
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
                if str(idx) in classifications:
                    for read, entry in classifications[str(idx)].items():
                        trusted_placements[base_file][read] = entry['placement']
        # Split the original jplace file
        # and write split jplaces to separate file directories
        with open(jplace) as f: jplace_json = json.load(f)
        cluster_dict = self.convert_cluster_dict_keys_to_aliases(clusterer.seq_library,
                                                                 alias_hash)
        hash_with_placements = self.jplace_split(jplace_json,
                                                 cluster_dict)

        for file_alias, placement_entries_list in hash_with_placements.items():
            alias_hash[file_alias]['place'] = placement_entries_list

        for k in alias_hash.keys():
            if 'place' not in alias_hash[k]:
                alias_hash[k]['place'] = []
        self.write_jplace(jplace_json,
                          alias_hash)

        self.hk.delete(files_to_delete)# Remove combined split, not really useful

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
