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

import sets

class TaxonomyCleaner:
    
    meaningless_taxonomic_names = sets.Set(['k__', 'd__', 'p__', 'c__', 'o__', 
                                            'f__', 'g__', 's__'])
    
    def remove_empty_ranks(self, tax_list):
        '''
        Removes empty rank prefixes
        
        Parameters
        ----------
        tax_list    : list
            A list of taxonomic ranks.
        Returns
        -------
        A list of taxonomic ranks with empty prefixes removed.
        '''
        new_tax_list = []
        for rank in tax_list:
            if rank not in self.meaningless_taxonomic_names:
                new_tax_list.append(rank)
        return new_tax_list
        