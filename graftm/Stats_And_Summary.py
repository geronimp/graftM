#!/usr/bin/env python


from Bio import SeqIO



class Stats_And_Summary:

    def check_read_length(self, reads, pipe):
        lengths = []
        
            
        for record in list(SeqIO.parse(open(reads, 'r'), 'fasta')):
            lengths.append(len(record.seq))
        if pipe == "P":
            return (sum(lengths) / float(len(lengths)))/3
        elif pipe =="D":
            return sum(lengths) / float(len(lengths))


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
                                    splt[2]
                                    ))
            
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

        
    def otu_builder(self, gup_file, output, cutoff, header):
        d = {}
        classifications = []
        placed = []
        otu_id = 0
        output_table = ['#ID\t'+header+'\tConsensusLineage']
        unique_list = []
    
    
        for line in open(gup_file, 'r'):
            list = line.split()
    
            if list[0] != 'name' and list[1] == list[2] and float(list[len(list)-2]) > float(cutoff):
    
                if list[0] not in d:
                    d[list[0]] = []
    
    
                d[list[0]].append(list[3])
    
        for x,y in d.iteritems():
    
            if x not in placed:
                classifications.append(';'.join(y))
                placed.append(x)
    
            else:
                continue
    
        for x in classifications:
    
            if x not in unique_list:
                unique_list.append(x)
    
        for x in unique_list:
            output_table.append([str(otu_id),str(classifications.count(x)),x])
            otu_id += 1
    
        with open(output, 'w') as otu_table:
    
            for line in output_table:
    
                if '#' in line:
                    otu_table.write(line+'\n')
    
                else:
                    otu_table.write('\t'.join(line)+'\n')

                    