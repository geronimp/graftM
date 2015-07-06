import numpy as np
from biom.table import Table
import tempfile
import subprocess
import logging

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

    def readTax(self, placements, output_path):
        with open(output_path, 'w') as out:
            for read, tax in placements.iteritems():
                out.write("%s\t%s\n" % (read, '; '.join(tax)))

    def build_basic_statistics(self, summary_hash, output, pipe):
        
        def compile_run_stats(hash):
            files = []
            tot_16S = []
            tot_18S = []
            cont_18S = []
            search_step = []
            placed_reads = 0
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
                    if hash['merge_reads']:
                        tot_16S += [str(len(hash[base]['forward']['reads'].keys()))]
                    else:
                        tot_16S += [str(len(hash[base]['reads'].keys()))]
                    try:
                        placed_reads+=int(tot_16S[0])-int(hash[base]['euk_contamination'])
                    except:

                        placed_reads+=int(tot_16S[0])
                    try:
                        tot_18S += [str(hash[base]['euk_contamination'] + hash[base]['euk_uniq'])]

                    except:
                        tot_18S += ['N/A']
                    try:
                        cont_18S += [str(summary_hash[base]['euk_contamination'])]
                    except:
                        cont_18S += ['N/A']
                    if hash['merge_reads']:
                        search_step += [hash[base]['forward']['search_t']]
                        aln_step += [hash[base]['forward']['aln_t']]
                        euk_check_step += [hash[base]['forward']['euk_check_t']]
                    else:
                        search_step += [hash[base]['search_t']]
                        aln_step += [hash[base]['aln_t']]
                        euk_check_step += [hash[base]['euk_check_t']]
                else:
                    raise Exception('Programming Error')
            return '\t'.join(files), '\t'.join(tot_16S), '\t'.join(tot_18S), '\t'.join(cont_18S), str(placed_reads), '\t'.join(files), '\t'.join(search_step), '\t'.join(aln_step), '\t'.join(euk_check_step), hash['place_t'], hash['summary_t'], hash['all_t']
        

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

    def _iterate_otu_table_rows(self, read_taxonomies):
        '''yield that which is required for an OTU table: taxonomy, and
        count of that taxonomy in each sample as an array
        
        Parameters
        ---------
        read_taxonomies:
            a list of hashes, where the position in the list corresponds to the 
            sample list, the key is the read name, and the value is an array
            of taxonomic info
            
        Yield
        -----
        OTU_ID: int, made up by this function
        taxonomy: array where each element is a taxonomic level
        count: number of observations of that taxonomy in each sample respectively
            
        Return
        ------
        Nothing, use this as an iterator'''
        
        taxonomy_string_to_taxonomy_array = {}
        taxonomy_string_to_counts = {}
        sample_index = 0
        num_samples = len(read_taxonomies)
        for read_to_taxonomy in read_taxonomies: # For each sample
            for taxonomy_array in read_to_taxonomy.values(): # For each read
                taxonomy_string = '; '.join(taxonomy_array)
                if taxonomy_string_to_taxonomy_array.has_key(taxonomy_string):
                    if taxonomy_string_to_taxonomy_array[taxonomy_string] != taxonomy_array:
                        raise Exception("Programming error: two different taxonomies had same taxonomy string")
                    try:
                        taxonomy_string_to_counts[taxonomy_string][sample_index] += 1
                    except KeyError:
                        taxonomy_string_to_counts[taxonomy_string][sample_index] = 1
                else:
                    taxonomy_string_to_taxonomy_array[taxonomy_string] = taxonomy_array
                    taxonomy_string_to_counts[taxonomy_string] = [0]*num_samples
                    taxonomy_string_to_counts[taxonomy_string][sample_index] = 1
            sample_index += 1
            
        otu_id = 1
        for tax_string, counts_array in taxonomy_string_to_counts.iteritems():
            array = []
            for i in range(num_samples):
                try:
                    array.append(counts_array[i])
                except IndexError:
                    array.append(0)
            yield otu_id,\
                taxonomy_string_to_taxonomy_array[tax_string],\
                array
            otu_id += 1
            
    def write_biom(self, sample_names, read_taxonomies, biom_file_io):
        '''Write the OTU info to a biom IO output stream
        
        Parameters
        ----------
        sample_names: String
            names of each sample (sample_ids for biom)
        read_taxonomies: Array of hashes as per _iterate_otu_table_rows()
        biom_file_io: io
            open writeable stream to write biom contents to
            
        Returns True if successful, else False'''
        counts = []
        observ_metadata = []
        otu_ids = []
        for otu_id, tax, count in self._iterate_otu_table_rows(read_taxonomies):
            if len(count) != len(sample_names):
                raise Exception("Programming error: mismatched sample names and counts")
            counts.append(count)
            observ_metadata.append({'taxonomy': tax})
            otu_ids.append(str(otu_id))
        
        table = Table(np.array(counts),
                      otu_ids, sample_names, observ_metadata,
                      [{}]*len(sample_names), table_id='GraftM Taxonomy Count Table')
        try:
            table.to_hdf5(biom_file_io, 'GraftM graft')
            return True
        except RuntimeError as e:
            logging.warn("Error writing BIOM output, file not written. The specific error was: %s" % e)
            return False

    def write_tabular_otu_table(self, sample_names, read_taxonomies, combined_output_otu_table_io):
        '''A function that takes a hash of trusted placements, and compiles them
        into an OTU-esque table.'''
        delim = u'\t'
        combined_output_otu_table_io.write(delim.join(['#ID', 
                                                       delim.join(sample_names),
                                                       'ConsensusLineage']))
        combined_output_otu_table_io.write(u"\n")
        for otu_id, tax, counts in self._iterate_otu_table_rows(read_taxonomies):
            combined_output_otu_table_io.write(delim.join(\
                (str(otu_id),
                 delim.join([str(c) for c in counts]),
                 '; '.join(tax)))+"\n")
            
    def write_krona_plot(self, sample_names, read_taxonomies, output_krona_filename):
        '''Creates krona plot at the given location. Assumes the krona executable
        ktImportText is available on the shell PATH'''
        tempfiles = []
        for n in sample_names:
            tempfiles.append(tempfile.NamedTemporaryFile(prefix='GraftMkronaInput', suffix=n))
        
        delim=u'\t'
        for _, tax, counts in self._iterate_otu_table_rows(read_taxonomies):
            for i, c in enumerate(counts):
                if c != 0:
                    tempfiles[i].write(delim.join((str(c),
                                                  delim.join(tax)
                                                  ))+"\n")
                    
        for t in tempfiles:
            t.flush()
        
        cmd = ["ktImportText",'-o',output_krona_filename]
        for i, tmp in enumerate(tempfiles):
            cmd.append(','.join([tmp.name,sample_names[i]]))

        # run the actual krona
        subprocess.check_output(' '.join(cmd), shell=True)

        # close tempfiles
        for t in tempfiles:
            t.close()

        
        