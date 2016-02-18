#!/usr/bin/env python
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
# System imports
import logging 

###############################################################################
###############################################################################
###############################################################################

class MalformedGreenGenesFormattedFile(Exception):
    pass

#-------------------------------------------------------------------------------

class ProgrammingError(Exception):
    pass

#-------------------------------------------------------------------------------

class MissingTaxonomy(Exception):
    pass

#-------------------------------------------------------------------------------

class Taxonomy:

    TAB='\t'
    COM=','
    GG_SEP_0='; '
    GG_SEP_1=';'
    
    def __init__(self):
        self.taxonomy = {}
        
    def _amend_taxonomy_gaps(self):
        '''
        Checks that each value in the self.taxonomy dict has consistent
        length. 
        '''
        self.max = max([len(x) for x in self.taxonomy.values()])
        for taxa_id, taxa_taxonomy_list in self.taxonomy.iteritems():
            try:
                if len(taxa_taxonomy_list) == self.max:
                    continue
                else:
                    ammended_split_taxonomy = \
                                taxa_taxonomy_list + \
                                [''] * (self.max - len(taxa_taxonomy_list))
                self.taxonomy[taxa_id] = ammended_split_taxonomy
            except:
                raise ProgrammingError("Error when trying to interpret taxonomy \
input on this line: %s\t%s" % (taxa_id, taxa_taxonomy_list))
    
    def __get__(self, instance, owner):
        '''
        Check whether an id is within the taxonomy description.
        
        Parameters
        ----------
        taxa_id : array
            taxa id to return the taxonomy for
        '''
        return self.taxonomy
        
    def entry_ids(self):
        '''
        Return a list of all taxonomy entries within the self.taxnomy dict
        '''
        return self.taxonomy.keys()
    
    def rank(self, taxa_id, rank):
        '''
        Parameters
        ----------
        taxa_id : array
            taxa id to return the taxonomy for
        rank : int
            rank to return
        Returns
        -------
        rank as indicated by the index passed to 'rank'.
        '''
        if rank < self.max:
            return self.taxonomy[taxa_id][rank]
        else:
            raise ProgrammingError("The rank specified was larger than the \
maximum possible number of ranks within the taxonomy description.")
        
        
    
    def taxonomy_string(self, taxa_id):
        '''
        Return the tax string of the id provided
        
        Parameters
        ----------
        taxa_id : array
            taxa id to return the taxonomy for
        
        Returns
        -------
        The taxonomy string joined with the standard self.GG_SEP_0 separator.
        '''
        
        if taxa_id in self.taxonomy:
            taxonomy_string = GG_SEP_0.join(self.taxonomy[taxa_id])
        else:
            raise MissingTaxonomy("Taxa id provided did not have any \
any corresponding entry in the taxonomy description.")
            
        return taxonomy_string
            
    def read_format_greengenes(self, taxonomy_io):
        '''
        Parses a greengenes taxonomy file. No set number of taxonomic ranks is 
        required
        
        Parameters
        ----------
        taxonomy_io : Open greengenes taxonomy file
        '''
        logging.debug("Reading in greengenes formatted taxonomy")
        for entry in taxonomy_io:
            try:                   
                if len(entry.strip()) == 0: continue #ignore empty lines 
                splits = entry.split(self.TAB)
                if len(splits) != 2:
                    raise MalformedGreenGenesTaxonomyException("Unexpected number of tab-separated fields found in taxonomy file, on line %s" % entry)
                name = splits[0].strip()
                taxonomy = [t.strip() for t in splits[1].split(self.GG_SEP_1)]
                while len(taxonomy) > 0 and not taxonomy[-1]:
                    taxonomy = taxonomy[:-1]
                for lineage in taxonomy:
                    if not lineage:
                        raise MalformedGreenGenesTaxonomyException("Encountered a taxonomy string with the middle of the taxonomy string missing: %s" % entry)
                
                if name in self.taxonomy:
                    raise MalformedGreenGenesTaxonomyException("Duplicate definition of taxonomy for %s" % name)
                self.taxonomy[name] = taxonomy
                
            except:
                raise MalformedGreenGenesFormattedFile("The GreenGenes taxonomy file \
provided was misformatted")
        self._amend_taxonomy_gaps()
        

            
    
    def read_format_taxtastic(self, taxonomy_io, seqinfo_io):
        logging.debug("Reading in taxtastic formatted taxonomy")
        raise ProgrammingError("Not implemented yet.")
      