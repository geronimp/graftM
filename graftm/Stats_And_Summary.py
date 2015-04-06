from Bio import SeqIO

class Stats_And_Summary:
    
    def __init__(self): pass

    def coverage_of_hmm(self, hmm, count_table, coverage_table, avg_read_length):
        
        for line in open(hmm):
            
            if line.startswith('LENG'):
                length = line.split()[1] 
        
        with open(coverage_table, 'w') as ct:
            write = ['#ID     20110816_S1D    ConsensusLineage\n']
            for line in open(count_table):
                    
                    if line.startswith('#'):
                        continue
                    
                    splt = line.split()
    
                    cov = str(round((float(splt[1])*float(avg_read_length)) / float(length), 3))
                        
                    write.append("%s\t%s\t%s\n" % (
                                    splt[0],
                                    cov,
                                    splt[2]))
            for entry in write:
                ct.write(entry)
            

    def build_basic_statistics(self, summary_hash, base_list, output, pipe):
        sep = '\t'
        S16= []
        S18 = []
        cont_18S = []
        n_placed = []
        search_s = []
        extract_t = []
        aln_t = []
        euk_check_t = []
        tree_i_t = summary_hash['place_t']
        summary_t = summary_hash['summary_t']
        all_t = str(round(summary_hash['stop_all'] - summary_hash['start_all'], 2))
        base_titles = '\t'.join(base_list)
        
        for base in base_list:
            h = summary_hash[base]
            if pipe == "P":
                S16.append(str(h['n_total_reads']))
                S18.append('N/A')
                cont_18S.append('N/A')
            elif pipe == "D":
                S16.append(str(h['n_total_reads'] - h['n_contamin_euks']))
                S18.append(str(h['n_uniq_euks']))
                cont_18S.append(str(h['n_contamin_euks']))
            
            n_placed.append(str(h['reads_found']))
            search_s.append(str(h['search_t']))
            extract_t.append(str(h['extract_t']))
            aln_t.append(str(h['extract_t']))
            euk_check_t.append(str(h['euk_check_t']))
        


        stats = """Basic run statistics (count):

                             Files:\t%s
Total number of 16S reads detected:\t%s
Total number of 18S reads detected:\t%s  
    'Contaminant' eukaryotic reads:\t%s
Number of 16S reads placed in tree:\t%s
        
Runtime (seconds):
                             Files:\t%s
                       Search step:\t%s
                      Extract step:\t%s
                    Alignment step:\t%s
              Eukaryote check step:\t%s

               Tree insertion step:\t%s
                 Summarising steps:\t%s
                     Total runtime:\t%s
    """ % (
                base_titles,
                sep.join(S16),
                sep.join(S18),
                sep.join(cont_18S),
                sep.join(n_placed),
                base_titles,
                sep.join(search_s),
                sep.join(extract_t),
                sep.join(aln_t),
                sep.join(euk_check_t),
                tree_i_t,
                summary_t,
                all_t
               )
        with open(output, 'w') as stats_file:
            stats_file.write(stats)

        
    def otu_builder(self, hash, output, base):
        ## A function that takes a hash of trusted placements, and compiles them
        ## into an OTU-esque table.
        # Start an output list and hash that will contain consolidated numbers 
        # for each taxonomic rank...
        output_table = [['#ID', base, 'ConsensusLineage']]
        c_hash = {}
        rank_id = 0
        # Consolidate the placement list into a simplified hash with taxonomic 
        # rank, and number. Simple.
        for read_name, placement_list in hash.iteritems():
            try:
                c_hash['; '.join(placement_list)] += 1
            except:
                c_hash['; '.join(placement_list)] = 1
        # Make an entry into the output list
        for rank, count in c_hash.iteritems():
            output_table.append([str(rank_id), str(count), rank])
            rank_id += 1
        # And write to a file
        with open(output, 'w') as out:
            for line in output_table:
                out.write('\t'.join(line)+'\n')