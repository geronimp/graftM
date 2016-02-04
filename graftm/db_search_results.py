#!/usr/bin/env python
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import logging
import re
import subprocess

from graftm.sequence_search_results import SequenceSearchResult
from _ctypes import alignment

class DBSearchResult:
    '''
    DBSearchResult - Class for containing results from search pipeline in GraftM
    '''
    HMM_NAME_FIELD = "NAME"
    SLASH_ENDING_REGEX = re.compile('.*/[12]$')
    ORFM_REGEX = re.compile('^(\S+)_(\d+)_(\d)_(\d+)')
    PREVIOUS_SPAN_CUTOFF = 0.25

    def __init__(self, gpkg_hash):
        '''
        Parameters
        ----------
        gpkg_hash : hash
            A hash containing the complete path to the graftm package as the 
            key and the GraftMPackage object as the entry
        '''
        
        self.alignment_hmms_hash = {}
        self.maximum_range_hash = {}
        self.hmm_name_hash = {}
        
        for gpkg in gpkg_hash.values():
            alignment_hmm = gpkg.alignment_hmm_path()
            search_hmms = gpkg.search_hmm_paths()
            for hmm in search_hmms:
                self.alignment_hmms_hash[hmm] = alignment_hmm
                for line in open(hmm):
                    if line.startswith(self.HMM_NAME_FIELD):
                        _, hmm_name = line.strip()\
                                          .split()     
                        self.hmm_name_hash[hmm_name] = hmm
                if gpkg.maximum_range():
                    self.maximum_range_hash[hmm] = gpkg.maximum_range()
                else:
                    raise Exception("NO MAX RANGE FOUND")
            
    def _get_hit_spans(self, search_result, maximum_range):
        '''
        _get_read_names - loops through hmm hits and their alignment spans to
        determine if they are potentially linked (for example, if one gene in a 
        contig hits a hmm more than once, in two different conserved regions
        of that gene) and combines them into one 'hit' because they are
        technically the same. The total span of the hits deemed to be linked is
        returned.

        Parameters
        ----------
        maximum_range : int
            Maximum range that a gene can extend within a contig. Any hits
            that extend beyond this length cannot be linked. max_range is
            set as 1.5 X the average length of all full length genes used
            in the search database. This is defined in the CONTENTS.json file
            within a gpkg.
        Returns
        -------
            Dictionary where keys are the contig/read name. The value for each
            entry is an array lists, one per hit in each contig, each with the
            span (min and max) of the alignment.
        '''

        splits = {}  # Define an output dictionary to be filled
        spans = list(
                     search_result.each(
                                        [SequenceSearchResult.QUERY_ID_FIELD,
                                         SequenceSearchResult.ALIGNMENT_DIRECTION,
                                         SequenceSearchResult.HIT_FROM_FIELD,
                                         SequenceSearchResult.HIT_TO_FIELD,
                                         SequenceSearchResult.QUERY_FROM_FIELD,
                                         SequenceSearchResult.QUERY_TO_FIELD]
                                        )
                    )

        for hit in spans:  # For each of these rows (i.e. hits)
            i = hit[0]  # set id to i
            c = hit[1]  # set complement to c
            ft = [min(hit[2:4]), max(hit[2:4])]  # set span as ft (i.e. from - to) - This is the amount covering the query
            qs = [min(hit[4:6]), max(hit[4:6])]  # seq the query span to qs - This is the amount covering the HMM

            if ft[0] == ft[1]: continue  # if the span covers none of the query, skip that entry (seen this before)

            if i not in splits:  # If the hit hasnt been seen yet
                splits[i] = {'span'       : [ft],
                             'strand'     : [c],
                             'query_span' : [qs]}  # add the span and complement as new entry

            else:  # otherwise (if it has been seen)
                for idx, entry in enumerate(splits[i]['span']):  # for each previously added entry
                    if splits[i]['strand'][idx] != c:  # If the hit is on the same complement strand
                        splits[i]['span'].append(ft)  # Add the new range to be split out in the future
                        splits[i]['strand'].append(c)  # Add the complement strand as well
                        splits[i]['query_span'].append(qs)
                        break

                    previous_qs      = splits[i]['query_span'][idx] # Get the query span of the previous hit

                    previous_q_range = set(range(previous_qs[0], previous_qs[1])) # Get the range of each
                    current_q_range  = set(range(qs[0], qs[1]))
                    query_overlap    = set(previous_q_range).intersection(current_q_range) # Find the intersection between the two ranges

                    previous_ft_span = set(range(entry[0], entry[1]))
                    current_ft_span  = set(range(ft[0], ft[1]))
                    if any(query_overlap): # If there is an overlap
                        ####################################################
                        # if the span over the actual read that hit the HMM
                        # for each hit overlap by > 25%, they are considered
                        # the same hit, and ignored
                        ####################################################
                        intersection_fraction = float(len(previous_ft_span.intersection(current_ft_span)))

                        if intersection_fraction / float(len(previous_ft_span)) >= self.PREVIOUS_SPAN_CUTOFF:
                            break
                        elif intersection_fraction / float(len(current_ft_span)) >= self.PREVIOUS_SPAN_CUTOFF:
                            break
                        else: # else (i.e. if the hit covers less that 25% of the sequence of the previous hit)
                            ####################################################
                            # If they made it this far, it means that the hits do not overlap.
                            # But one last check must be made to ensure they do not cover the same
                            # region in the HMM.
                            ####################################################
                            if len(query_overlap) > (len(current_q_range)*PREVIOUS_SPAN_CUTOFF): # if the overlap on the query HMM does not span over 25%
                                if (idx+1) == len(splits[i]['span']):
                                    splits[i]['span'].append(ft) # Add from-to as another entry, this is another hit.
                                    splits[i]['strand'].append(c) # Add strand info as well
                                    splits[i]['query_span'].append(qs)
                                    break # And break

                    if min(entry) < min(ft):  # if/else to determine which entry comes first (e.g. 1-5, 6-10 not 6-10, 1-5)
                        if max(ft) - min(entry) < maximum_range:  # Check if they lie within range of eachother
                            entry[1] = max(ft)  # ammend the entry if they are
                            break  # And break the loop
                    else:
                        if max(entry) - min(ft) < maximum_range:  # Check if they lie within range of eachother
                            entry[0] = min(ft)  # ammend the entry if they are
                            break  # And break the loop

                else:  # if no break occured (no overlap)
                    splits[i]['span'].append(ft)  # Add the new range to be split out in the future
                    splits[i]['strand'].append(c)  # Add the complement strand as well
                    splits[i]['query_span'].append(qs)
        
        return {key: {"entry":entry['span'], 'strand': entry['strand']} for key, entry in splits.iteritems()}

    def read(self, search_results):
        '''
        Parameters
        ----------
        Read in the search results object from SequenceSearchResult objects
        
        search_results: obj
            SequenceSearchResult object from HMM or DIAMOND pipeline.
        '''
        logging.debug("Parsing search results")        
        self.hits={}
        self.alignment_batches = {}
        for search_result in search_results:
            self.search_hmm_name=None
            for read_data in search_result.each([search_result.QUERY_ID_FIELD,
                                                 search_result.HMM_NAME_FIELD]): 
                hit_object = HitObject(read_data)
                if hit_object.clean_name() in self.hits:
                    raise Exception("Non-redundant read names were provided \
in the sequence files. Please give each read a unique name")
                else:
                    self.hits[hit_object.name()] = hit_object
                if self.search_hmm_name:
                    if read_data[1]!=self.search_hmm_name:
                        raise Exception("Error parsing results: \
hmmsearch or diamond search results have multiple entries for search hmm.")   
                else:
                    self.search_hmm_name = read_data[1]
                alignment_hmm = self.alignment_hmms_hash[
                                    self.hmm_name_hash[self.search_hmm_name]]
                if alignment_hmm in self.alignment_batches:
                    self.alignment_batches[alignment_hmm]\
                                                .append(hit_object.name())
                else:
                    self.alignment_batches[alignment_hmm] = [hit_object.name()]

            maximum_range = \
                self.maximum_range_hash[self.hmm_name_hash[self.search_hmm_name]]
            spans = self._get_hit_spans(search_result, 
                                        maximum_range)        
            self.spans
    
class HitObject(DBSearchResult):
    
    def __init__(self, hit):
        self.read_name_string, \
        self.hmm_name_string = hit
        self.clean_name_string = self.read_name_string
        self._find_orfm_suffix()
        self._find_slash()

    def _find_orfm_suffix(self):
        self.orfmregex = \
            self.ORFM_REGEX.match(self.read_name_string)
        if self.orfmregex:
            self.clean_name_string = self.orfmregex.groups(0)[0] 
            
    def _find_slash(self):
        self.slashendingregex = \
            self.SLASH_ENDING_REGEX.match(self.clean_name_string)
        if self.slashendingregex:
            self.clean_name_string = self.clean_name_string[:-2]
            
    def clean_name(self):
        return self.clean_name_string
    
    def name(self):
        return self.read_name_string
