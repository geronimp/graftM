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
################################# - Imports - ##################################

# Local imports
from graftm.getaxnseq import Getaxnseq
from graftm.greengenes_taxonomy import GreenGenesTaxonomy, MalformedGreenGenesTaxonomyException
from graftm.taxonomy_cleaner import TaxonomyCleaner

# System imports
import logging

################################################################################
################################## - Code - ####################################

class TreeDecorator:
    '''
    A class that conservatively decorates trees with taxonomy, or any other
    hierarchical annotation. If all tips descending from a node within the
    provided tree have consistent taxonomy, it will be decorated with that
    taxonomy (or annotation of any type).
    '''

    def __init__(self, tree, taxonomy, seqinfo=None):
        '''
        Parameters
        ----------
        tree        : dendropy.Tree

            dendropy.Tree object
        taxonomy    : string
            Path to a file containing taxonomy information about the tree,
            either in Greengenes or taxtastic format (seqinfo file must also
            be provided if taxonomy is in taxtastic format).
        seqinfo     : string
            Path to a seqinfo file. This is a .csv file with the first column
            denoting the sequence name, and the second column, its most resolved
            taxonomic rank.
        '''

        self.encountered_nodes = {}
        self.encountered_taxonomies = set()
        self.tree = tree

        # Read in taxonomy
        logging.info("Reading in taxonomy")
        if seqinfo:
            logging.info("Importing taxtastic taxonomy from files: %s and %s" % (taxonomy, seqinfo))
            gtns = Getaxnseq()
            self.taxonomy =  gtns.read_taxtastic_taxonomy_and_seqinfo(open(taxonomy), open(seqinfo))

        else:
            try:
                logging.info("Reading Greengenes style taxonomy")
                self.taxonomy = GreenGenesTaxonomy.read_file(taxonomy).taxonomy
            except MalformedGreenGenesTaxonomyException:
                raise Exception("Failed to read taxonomy as a Greengenes \
                                 formatted file. Was a taxtastic style \
                                 taxonomy provided with no seqinfo file?")

    def _write_consensus_strings(self, output):
        '''
        Writes the taxonomy of each leaf to a file. If the leaf has no
        taxonomy, a taxonomy string will be created using the annotations
        provided to the ancestor nodes of that leaf (meaning, it will be
        decorated).

        Parameters
        ----------
        output    : string
            File to which the taxonomy strings for each leaf in the tree will
            be written in Greengenes format, e.g.
                637960147    mcrA; Euryarchaeota_mcrA; Methanomicrobia
                637699780    mcrA; Euryarchaeota_mcrA; Methanomicrobia
        '''

        logging.info("Writing decorated taxonomy to file: %s" % (output))

        with open(output, 'w') as out:

            for tip in self.tree.leaf_nodes():
                tax_name = tip.taxon.label.replace(" ", "_")

                if tip.taxon.label in self.taxonomy:
                    tax_string = '; '.join(self.taxonomy[tax_name])
                else:
                    ancestor_list = []
                    for ancestor in tip.ancestor_iter():
                        if ancestor.label:
                            split_node_name = ancestor.label.split(':')
                            if len(split_node_name) == 2:
                                ancestor_list+=list(reversed(split_node_name[1].split('; ')))
                            elif len(split_node_name) == 1:
                                try:
                                    float(split_node_name[0])
                                except ValueError:
                                    ancestor_list+=list(reversed(split_node_name[0].split('; ')))
                            else:
                                raise Exception("Malformed node name: %s" % ancestor.label)
                    tax_list = list(reversed(ancestor_list))

                    if len(tax_list) < 1:
                        logging.warning("No taxonomy found for species %s!" % (tax_name))
                        tax_string = "Unknown"
                    else:
                        tax_string = '; '.join(tax_list)

                output_line = "%s\t%s\n" % (tax_name, tax_string)
                out.write(output_line)

    def _rename(self, node, name):
        '''
        Rename an internal node of the tree. If an annotation is already
        present, append the new annotation to the end of it. If a bootstrap
        value is present, add annotations are added after a ":" as per standard
        newick format.

        Parameters
        ----------
        node: dendropy.Node
            dendropy.Node object
        name    : string
            Annotation to rename the node with.
        '''
        if node.label:
            try:
                float(node.label)
                new_label = "%s:%s" % (node.label,
                                       name)
            except ValueError:
                new_label = "%s; %s" % (node.label,
                                        name)
            node.label = new_label
        else:
            node.label = name


    def decorate(self, output_tree, output_tax, unique_names):
        '''
        Decorate a tree with taxonomy. This code does not allow inconsistent
        taxonomy within a clade. If one sequence in a clade has a different
        annotation to the rest, it will split the clade. Paraphyletic group
        names are distinguished if unique_names = True using a simple tally of
        each group (see unique_names below).

        Parameters
        ----------
        output_tree        : string
            File to which the decorated tree will be written.
        output_tax         : string
            File to which the taxonomy strings for each tip in the tree will be
            written.
        unique_names       : boolean
            True indicating that a unique number will be appended to the end of
            a taxonomic rank if it is found more than once in the tree
            (i.e. it is paraphyletic in the tree). If false, multiple clades
            may be assigned with the same name.
        '''
        logging.info("Decorating tree")
        encountered_taxonomies = {}
        tc = TaxonomyCleaner()
        for node in self.tree.preorder_internal_node_iter(exclude_seed_node=True):

            max_tax_string_length = 0
            for tip in node.leaf_nodes():
                tip_label=tip.taxon.label.replace(' ', '_')
                if tip_label in self.taxonomy:
                    tax_string_length \
                        = len(self.taxonomy[tip.taxon.label.replace(' ', '_')])
                    if tax_string_length > max_tax_string_length:
                        max_tax_string_length = tax_string_length

            logging.debug("Number of ranks found for node: %i" % max_tax_string_length)
            tax_string_array = []
            for rank in range(max_tax_string_length):

                rank_tax = []
                for tip in node.leaf_nodes():

                    tip_label = tip.taxon.label.replace(' ', '_')
                    if tip_label in self.taxonomy:

                        tip_tax = self.taxonomy[tip_label]
                        if len(tip_tax) > rank:
                            tip_rank = tip_tax[rank]
                            if tip_rank not in rank_tax:
                                rank_tax.append(tip_rank)
                consistent_taxonomy = len(rank_tax) == 1

                if consistent_taxonomy:
                    tax=rank_tax.pop()
                    logging.debug("Consistent taxonomy found for node: %s" \
                                                                % tax)
                    if tax not in tc.meaningless_taxonomic_names:
                        if unique_names:
                            if tax in encountered_taxonomies:
                                encountered_taxonomies[tax]+=0
                                tax = "%s_%i" \
                                        % (tax, encountered_taxonomies[tax])
                            else:
                                encountered_taxonomies[tax]=0
                        tax_string_array.append(tax)

            if any(tax_string_array):
                index = 0
                for anc in node.ancestor_iter():
                    try:
                        index+=anc.tax
                    except:
                        continue
                tax_string_array = tax_string_array[index:]
                if any(tax_string_array):
                    self._rename(node, '; '.join(tax_string_array))
                node.tax = len(tax_string_array)
        logging.info("Writing decorated tree to file: %s" % output_tree)
        if output_tree:
            self.tree.write(path=output_tree, schema="newick")
        if output_tax:
            self._write_consensus_strings(output_tax)

################################################################################
################################################################################
