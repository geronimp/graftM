import logging
from graftm.sequence_search_results import SequenceSearchResult

class SearchTableWriter:
    '''
    Class for writing the search output OTU table. Basically a summary
    of hits to the HMM/Diamond searched in the following format:
    
             #ID    Metagenome_1    Metagenome_2    ...
            HMM1    50              6
            HMM2    195             41
            HMM3    2               20120
            ...
            
    You just need to specify a series of SequenceSearchResult objects, and an
    output path.
    '''
    
    def _interpret_hits(self, results_list, base_list):
        '''Sort reads that hit multiple HMMs to the databases to which they had
        the highest bit score. Return a dictionary containing HMMs as keys, and 
        number of hits as the values.
        
        This function is set up so that the read names could easily be returned
        instead of numbers, for future development of GraftM
        
        Parameters
        ----------
        results_list: list
            Iterable if SequenceSearchResult objects. e.g.
                [SequenceSearchResult_1, SequenceSearchResult_2, ...]
        
        base_list: list
            Iterable of the basenames for each sequence file provided to graftM
            e.g.
                [sample_1, sample_2, ...]
                
        Returns
        -------
        dictionary:
            Contains samples as entries. The value for each sample is another 
            dictionary with HMM as the key, and number of hits as values:
                {"sample_1":{HMM_1: 12
                             HMM_2: 35
                             HMM_3: 1258
                             ...}
                 ...
                }
            
        '''
        logging.debug("Sorting reads into HMMs by bit score")
        
        run_results = {}

        ########################################################################
        ################## - Sort reads to best hit db - #######################
        for base, results in zip(base_list, results_list): # For each sample
            search_results = {}
            for search in results():
                search_list = list(
                                   search.each([SequenceSearchResult.QUERY_ID_FIELD,
                                                SequenceSearchResult.ALIGNMENT_BIT_SCORE,
                                                SequenceSearchResult.HMM_NAME_FIELD])
                                   )
                for hit in search_list:
                    if hit[0] in search_results:
                        if float(hit[1]) > search_results[hit[0]][1]:
                            search_results[hit[0]] = [float(hit[1]), hit[2]]
                    else:
                        search_results[hit[0]] = [float(hit[1]), hit[2]]
            run_results[base] = search_results
        
        ########################################################################
        ################## - Gather counts for each db - #######################
        db_count = {}
        for run in run_results.keys():
            run_count = {}
            for entry in run_results[run].values():
                key = entry[1]
                if key in run_count:
                    run_count[key] += 1
                else:
                    run_count[key] = 1
            db_count[run] = run_count
        
        return db_count    
    
    def _write_results(self, db_count, output_path):
        '''Write the table to the output_path directory
        
        db_count: dict
            Contains samples as entries. The value for each sample is another 
            dictionary with HMM as the key, and number of hits as values:
                {"sample_1":{HMM_1: 12
                             HMM_2: 35
                             HMM_3: 1258
                             ...}
                 ...
                }
        
        output_path: str
            Path to output file to which the resultant output file will be 
            written to.
        '''
        
        logging.debug("Writing search otu table to file: %s" % output_path)
        
        output_dict = {}
        
        for idx, value_dict in enumerate(db_count.values()):
            for database, count in value_dict.iteritems():
                if database in output_dict:
                    output_dict[database].append(str(count))
                else:
                    output_dict[database] = ['0']*idx + [str(count)]
            
            for key, item in output_dict.iteritems():
                if len(item) == idx:
                    output_dict[key].append('0')
        
        with open(output_path, 'w') as out:
            out.write('\t'.join(["#ID"] + db_count.keys()) + '\n')
            for key, item in output_dict.iteritems():
                out.write("%s\t%s" % (key, '\t'.join(item)) + '\n' )
        
    def build_search_otu_table(self, search_results_list, base_list, output_path):
        '''
        Build an OTU from SequenceSearchResult objects
        
        Parameters
        ----------
        search_results_list: list
            Iterable if SequenceSearchResult objects. e.g.
                [SequenceSearchResult_1, SequenceSearchResult_2, ...]
        base_list: list
            Iterable of the basenames for each sequence file provided to graftM
            e.g.
                [sample_1, sample_2, ...]
        output_path: str
            Path to output file to which the resultant output file will be 
            written to.
        '''
        
        db_count = self._interpret_hits(search_results_list,
                                        base_list)
        
        self._write_results(db_count, output_path)
        
        
        
