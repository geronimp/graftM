#!/usr/bin/env python

import re
import sets

import pprint as pp

def main(base, taxonomy_file):
    seqinfo = base+"_seqinfo.csv"
    output_tax = base+"_taxonomy.csv"
    
    Getaxnseq().gg_taxonomy_builder(taxonomy_file, output_tax, seqinfo)
    
    return seqinfo, output_tax

class Getaxnseq:

    def __init__(self): pass
    
    def taxonomy_line(self, level_index, taxon_array):
        if level_index == 0:
            return '%s,Root,kingdom,%s,%s,%s,,,,,,' % (taxon_array[level_index], taxon_array[level_index], 'Root', taxon_array[level_index])
        elif level_index == 1:
            return '%s,%s,phylum,%s,%s,%s,%s,,,,,' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1])
        elif level_index == 2:
            return '%s,%s,class,%s,%s,%s,%s,%s,,,,' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1], taxon_array[2])
        elif level_index == 3:
            return '%s,%s,order,%s,%s,%s,%s,%s,%s,,,' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1], taxon_array[2], taxon_array[3])
        elif level_index == 4:
            return '%s,%s,family,%s,%s,%s,%s,%s,%s,%s,,' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1], taxon_array[2], taxon_array[3], taxon_array[4])
        elif level_index == 5:
            return '%s,%s,genus,%s,%s,%s,%s,%s,%s,%s,%s,' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1], taxon_array[2], taxon_array[3], taxon_array[4], taxon_array[5])
        elif level_index == 6:
            return '%s,%s,species,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (taxon_array[level_index], taxon_array[level_index-1], taxon_array[level_index], 'Root', taxon_array[0], taxon_array[1], taxon_array[2], taxon_array[3], taxon_array[4], taxon_array[5], taxon_array[6])
        else:
            raise Exception("Programming error, found too many levels!")


    def gg_taxonomy_builder(self, taxonomy_file, output_taxonomy, output_seqinfo):
        first_pass_id_and_taxonomies = []
        meaningless_taxonomic_names = sets.Set(['k__', 'd__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'])

        for entry in open(taxonomy_file):
            # split the entire line including the  on '; '
            split = entry.rstrip().split('; ')

            # split out the taxon ID from the first split of above
            taxon_id, first_taxon = split[0].split()[:2]
            tax_split = [first_taxon] + split[1:]

            # Replace spaces with underscores e.g. 'Candidatus my_genus'
            for idx, item in enumerate(tax_split):
                tax_split[idx] = re.sub('\s+', '_', item)

            # Remove 'empty' taxononomies e.g. 's__'
            tax_split = [item for item in tax_split if item not in meaningless_taxonomic_names]

            # Add this fixed up list to the list
            first_pass_id_and_taxonomies.append([taxon_id]+tax_split)

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

        #pp.pprint(parents)
        #pp.pprint(first_pass_id_and_taxonomies)

        # Write the sequence file
        with open(output_seqinfo, 'w') as seqout:
            # write header
            seqout.write('seqname,tax_id\n')
            # write each taxonomic association
            for array in first_pass_id_and_taxonomies:
                seqout.write("%s,%s\n" % (array[0], array[-1]))

        # Write the taxonomy file
        noted_taxonomies = sets.Set([])
        with open(output_taxonomy, 'w') as seqout:
            # write header and root line
            seqout.write('tax_id,parent_id,rank,tax_name,root,kingdom,phylum,class,order,family,genus,species\n')
            seqout.write('Root,Root,root,Root,Root,,,,,,,\n')
            # write all the taxonomies
            for array in first_pass_id_and_taxonomies:
                taxons = array[1:]



                for i, tax in enumerate(taxons):
                    line = self.taxonomy_line(i, taxons[:(i+1)])

                    if line not in noted_taxonomies:
                        seqout.write(line+"\n")
                        noted_taxonomies.add(line)


