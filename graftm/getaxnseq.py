#!/usr/bin/env python
################################################################################
#                                                                              #
#     This program is free software: you can redistribute it and/or modify     #
#     it under the terms of the GNU General Public License as published by     #
#     the Free Software Foundation, either version 3 of the License, or        #
#     (at your option) any later version.                                      #
#                                                                              #
#     This program is distributed in the hope that it will be useful,          #
#     but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#     GNU General Public License for more details.                             #
#                                                                              #
#     You should have received a copy of the GNU General Public License        #
#     along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                              #
################################################################################

import re
import sets
import string
import logging

from graftm.taxonomy_cleaner import TaxonomyCleaner

class Getaxnseq:
    
    def _taxonomy_line(self, level_index, taxon_array, taxonomic_level_names):   
        if level_index == 0:
            return '%s,Root,%s,%s,%s,%s%s' \
                % (taxon_array[level_index],
                    #parent
                   taxonomic_level_names[level_index],
                   taxon_array[level_index],
                   'Root',
                   taxon_array[level_index],
                   ''.join([','] * (len(taxonomic_level_names) - (level_index+1)) ))
        else:
            return '%s,%s,%s,%s,%s,%s%s' \
                % (taxon_array[level_index],
                   taxon_array[level_index-1],
                   taxonomic_level_names[level_index],
                   taxon_array[level_index],
                   'Root',
                   ','.join(taxon_array),
                   ''.join([','] * (len(taxonomic_level_names) - (level_index+1)) ))
    
    def read_taxtastic_taxonomy_and_seqinfo(self, taxonomy_io, seqinfo_io):
        '''Read the taxonomy and seqinfo files into a dictionary of 
        sequence_name => taxonomy, where the taxonomy is an array of lineages
        given to that sequence.
        
        Possibly this method is unable to handle the full definition of these
        files? It doesn't return what each of the ranks are, for starters.
        Doesn't deal with duplicate taxon names either.
        '''
        
        # read in taxonomy file
        lineages = [] #array of hashes where each taxon points to its parent taxon's name
        taxon_to_lineage_index = {}
        expected_number_of_fields = None
        for line in taxonomy_io:
            splits = line.strip().split(',')
            if expected_number_of_fields is None:
                expected_number_of_fields = len(splits)
                lineages = [{}]* (expected_number_of_fields-4)
                continue #this is the header line
            elif len(splits) != expected_number_of_fields:
                raise Exception("Encountered error parsing taxonomy file, expected %i fields but found %i on line: %s" %
                                (expected_number_of_fields, len(splits), line))
            # e.g. 'tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species
            tax_id = splits[0]
            parent_id = splits[1]
            
            try:
                lineage_index = splits.index('')-5
            except ValueError:
                lineage_index = len(splits)-5
            taxon_to_lineage_index[tax_id] = lineage_index
            lineages[lineage_index][tax_id] = parent_id
        
        taxonomy_dictionary = {}
        for i, line in enumerate(seqinfo_io):
            if i==0: continue #skip header line
            
            splits = line.strip().split(',')
            if len(splits) != 2:
                raise Exception("Bad formatting of seqinfo file on this line: %s" % line)
            
            seq_name = splits[0]
            taxon = splits[1]
            lineage_index = taxon_to_lineage_index[taxon]
            if lineage_index==0:
                # Don't include Root in the taxonomy
                taxonomy_dictionary[seq_name] = []
            else:
                full_taxonomy_rev = []
                while lineage_index > 0:
                    full_taxonomy_rev.append(taxon)
                    taxon = lineages[lineage_index][taxon]
                    lineage_index = lineage_index-1
                taxonomy_dictionary[seq_name] = list(reversed(full_taxonomy_rev))
            
        return taxonomy_dictionary
    
    def write_taxonomy_and_seqinfo_files(self, taxonomies, output_taxonomy_file, output_seqinfo_file):
        '''Write out taxonomy and seqinfo files as required by taxtastic
        from known taxonomies
        
        Parameters
        ----------
        taxonomies:
            hash of taxon_id to array of taxonomic information
        output_taxonomy_file:
            write taxtastic-compatible 'taxonomy' file here
        output_seqinfo_file:
            write taxtastic-compatible 'seqinfo' file here'''
        
        first_pass_id_and_taxonomies = []
        tc=TaxonomyCleaner()
        max_number_of_ranks = 0
        
        for taxon_id, tax_split in taxonomies.iteritems():
            # Replace spaces with underscores e.g. 'Candidatus my_genus'
            for idx, item in enumerate(tax_split):
                tax_split[idx] = re.sub('\s+', '_', item.strip())

            # Remove 'empty' taxononomies e.g. 's__'
            tax_split = tc.remove_empty_ranks(tax_split)

            # Add this fixed up list to the list
            first_pass_id_and_taxonomies.append([taxon_id]+tax_split)

            if len(tax_split) > max_number_of_ranks:
                max_number_of_ranks = len(tax_split)

        # Find taxons that have multiple parents, building a hash of parents as we go (i.e. a tree of taxonomies embedded in a hash)
        #
        # Assumes that no two taxonomic labels are the same when they are from different
        # taxonomic levels. When there are children with multiple parents at the
        # same taxonomic label then these are warned about and worked around.
        parents = {} #hash of taxon to immediate parent
        known_duplicates = sets.Set([])
        for j, array in enumerate(first_pass_id_and_taxonomies):
            taxonomy = array[1:]
            for i, tax in enumerate(taxonomy):
                if i==0: continue #top levels don't have parents
                ancestry = taxonomy[i-1]
                if parents.has_key(tax):
                    if parents[tax] != ancestry:

                        dup = "%s%s" %(parents[tax], ancestry)
                        # don't report the same problem several times
                        if dup not in known_duplicates:
                            print " %s '%s' with multiple parents %s and %s" % (array[0], tax, parents[tax], ancestry)
                            known_duplicates.add(dup)

                        # fix the current one
                        new_name_id = 1
                        new_name = "%se%s" % (tax, new_name_id)

                        while parents.has_key(new_name) and parents[new_name] != ancestry:
                            new_name_id += 1
                            new_name = "%se%s" % (tax, new_name_id)

                        first_pass_id_and_taxonomies[j][i+1] = new_name
                        taxonomy[i] = new_name
                        parents[new_name] = ancestry
                else:
                    # normal case, seeing a new taxon and parent for the first time

                    parents[tax] = ancestry

        # Write the sequence file
        with open(output_seqinfo_file, 'w') as seqout:
            # write header
            seqout.write('seqname,tax_id\n')
            # write each taxonomic association
            for array in first_pass_id_and_taxonomies:
                sequence_id = array[0]
                if len(array)==1:
                    most_specific_taxonomic_affiliation = 'Root'
                else:  
                    most_specific_taxonomic_affiliation = array[-1]

                seqout.write("%s,%s\n" % (array[0], 
                                          most_specific_taxonomic_affiliation))

        # Write the taxonomy file
        noted_taxonomies = sets.Set([])
        taxonomic_level_names = ["rank_%i" % rank for rank in range(max_number_of_ranks)]
        with open(output_taxonomy_file, 'w') as seqout:
            # write header and root line
            seqout.write(','.join(['tax_id','parent_id','rank','tax_name','root'] + taxonomic_level_names) +'\n')
            seqout.write(','.join(['Root','Root','root','Root','Root'])+ ''.join([',']*max_number_of_ranks) + '\n')
            # write all the taxonomies
            for array in first_pass_id_and_taxonomies:
                taxons = array[1:]
                for i, tax in enumerate(taxons):
                    line = self._taxonomy_line(i, taxons[:(i+1)], taxonomic_level_names)
                    if line not in noted_taxonomies:
                        seqout.write(line+"\n")
                        noted_taxonomies.add(line)



