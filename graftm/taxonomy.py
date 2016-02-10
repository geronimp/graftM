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
__author__ = "Joel Boyd, Ben Woodcroft"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
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
        
    def _ammend_taxonomy_gaps(self):
        '''
        Checks that each value in the self.taxonomy dict has consistent
        length. 
        '''
        try:
            self.max = max([len(x) for x in self.taxonomy.values()])
            for taxa_id, taxa_taxonomy_list in self.taxonomy.iteritems():
                
                if len(taxa_taxonomy_list) == self.max:
                    continue
                else:
                    ammended_split_taxonomy = \
                                taxa_taxonomy_list + \
                                [''] * (self.max - len(taxa_taxonomy_list))
                self.taxonomy[taxa_id] = ammended_split_taxonomy
        except:
            logging.error("Controlled exit from taxonomy.Taxonomy \
 _ammend_taxonomy_gaps")    
            ProgrammingError("Error when trying to interpret taxonomy input")
    
    def check(self, taxa_id):
        '''
        Check whether an id is within the taxonomy description.
        
        Parameters
        ----------
        taxa_id : array
            taxa id to return the taxonomy for
        '''
        if taxa_id in self.id_set:
            return True
        else:
            return False
        
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
        try:
            for entry in taxonomy_io:
                entry_id, entry_taxonomy_string = entry.strip()\
                                                       .split(self.TAB)
                
                if self.GG_SEP_0 in entry_taxonomy_string:
                    split_entry_taxonomy = \
                                entry_taxonomy_string.split(self.GG_SEP_0)                
                else:
                    split_entry_taxonomy = \
                                 entry_taxonomy_string.split(self.GG_SEP_1)
                
                self.taxonomy[entry_id] = split_entry_taxonomy
            self._ammend_taxonomy_gaps()
            self.id_set = set(self.taxonomy.keys())
        except:
            logging.error("Controlled exit from taxonomy.Taxonomy \
read_format_greengenes")
            raise MalformedGreenGenesFormattedFile("The GreenGenes taxonomy file \
provided was misformatted")
    
    def read_format_taxtastic(self, taxonomy_io, seqinfo_io):
        logging.debug("Reading in taxtastic formatted taxonomy")
        raise ProgrammingError("Not implemented yet.")
      