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


    def build_basic_statistics(self, summary_hash, output, pipe):
        
        def compile_run_stats(hash):
            files = []
            tot_16S = []
            tot_18S = []
            cont_18S = []
            search_step = []
            placed_reads = []
            aln_step = []
            euk_check_step = []
            for base in hash['base_list']:
                if hash['reverse_pipe']:
                    files += [base + ' forward', base + ' reverse']
                    tot_16S += [str(len(hash[base]['forward']['reads'].keys())), str(len(hash[base]['reverse']['reads'].keys()))]
                    try:
                        tot_18S += [str(hash[base]['forward']['euk_uniq'] + hash[base]['forward']['euk_contamination']), 
                                    str(hash[base]['reverse']['euk_uniq'] + hash[base]['reverse']['euk_contamination'])]
                    except:
                        tot_18S += ['N/A', 'N/A']
                    try:
                        cont_18S += [str(summary_hash[base]['forward']['euk_contamination']), str(summary_hash[base]['reverse']['euk_contamination'])] 
                    except:
                        cont_18S = ['N/A']
                    search_step += [hash[base]['forward']['search_t'], hash[base]['reverse']['search_t']]
                    aln_step += [hash[base]['forward']['aln_t'], hash[base]['reverse']['aln_t']]
                    euk_check_step += [hash[base]['forward']['euk_check_t'], hash[base]['reverse']['euk_check_t']]
                elif not hash['reverse_pipe']:
                    files += [base]
                    tot_16S += [str(len(hash[base]['reads'].keys()))]
                    try:
                        placed_reads=[str(int(tot_16S[0])-int(hash[base]['euk_contamination']))]
                    except:
                        placed_reads=tot_16S
                    try:
                        tot_18S += [str(hash[base]['euk_contamination'] + hash[base]['euk_uniq'])]
                    except:
                        tot_18S += ['N/A']
                    try:
                        cont_18S += [str(summary_hash[base]['euk_contamination'])]
                    except:
                        cont_18S += ['N/A']
                    search_step += [hash[base]['search_t']]
                    aln_step += [hash[base]['aln_t']]
                    euk_check_step += [hash[base]['euk_check_t']]
                else:
                    raise Exception('Programming Error')
            return '\t'.join(files), '\t'.join(tot_16S), '\t'.join(tot_18S), '\t'.join(cont_18S), '\t'.join(placed_reads), '\t'.join(files), '\t'.join(search_step), '\t'.join(aln_step), '\t'.join(euk_check_step), hash['place_t'], hash['summary_t'], hash['all_t']
        
        stats = """Basic run statistics (count):

                             Files:\t%s
Total number of 16S reads detected:\t%s
Total number of 18S reads detected:\t%s
    'Contaminant' eukaryotic reads:\t%s
              reads placed in tree:\t%s

Runtime (seconds):
                             Files:\t%s
                       Search step:\t%s
                    Alignment step:\t%s
              Eukaryote check step:\t%s

               Tree insertion step:\t%s
                 Summarising steps:\t%s
                     Total runtime:\t%s
    """ % compile_run_stats(summary_hash)
        
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

