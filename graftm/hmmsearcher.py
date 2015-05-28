import subprocess
import logging

class HmmSearcher:
    r"""Runs hmmsearch given one or many HMMs in a scalable and fast way""" 
    
    def __init__(self, num_cpus, extra_args=''):
        r"""New
        
        Parameters
        ----------
        num_cpus: Integer
            The total number of CPUs to use when searching
        extra_args: String
            Extra arguments for hmmsearch. --cpus is already defined, so
            do not specify it in this argument"""
        self._num_cpus = num_cpus
        self._extra_args = extra_args
        
    def hmmsearch(self, input_pipe, hmms, output_files):
        r"""Run HMMsearch with all the HMMs, generating output files
        
        Parameters
        ----------
        input_pipe: String
            A string which is a partial command line. When this command is run
            is outputs to STDOUT fasta formatted protein sequences, which 
            hmmsearch runs on.
        hmms: list of paths
            A list of (string) paths to HMM files which are used to search with.
        output_files: list of paths
            A list of (string) paths to output CSV files to be generated by the
            HMM searching
            
        Returns
        -------
        N/A
        
        May raise an exception if hmmsearching went amiss"""
        
        # Check input and output paths are the same length
        if len(hmms) != len(output_files):
            raise Exception("Programming error: number of supplied HMMs differs from the number of supplied output files")
        
        # Create queue data structure
        queue = []
        for i, hmm in enumerate(hmms):
            queue.append( [hmm, output_files[i]] )
        
        # While there are more things left in the queue
        while len(queue) > 0:
            pairs_to_run = self.__munch_off_batch(queue)
            
            # Run hmmsearches with each of the pairs
            cmd = self.__hmm_command(input_pipe, pairs_to_run)
            logging.debug("Running command: %s" % cmd)
            subprocess.check_call(['bash','-c', cmd])
            
    def __munch_off_batch(self, queue):
        r"""Take a batch of sequences off the queue, and return pairs_to_run.
        The queue given as a parameter is affected
        """
        
        # if the number of CPUs used == 1, just pop one off (which == below)
        # elif the the number of things in the queue is 1, use all the CPUs and just that hmm
        if len(queue) == 1 or self._num_cpus == 1:
            pairs_to_run = [[queue.pop(0), self._num_cpus]]
        else:
            # else use 1 CPU per HMM.
            # TODO: if len(queue) < num_cpus, use more CPUs from some HMMs
            pairs_to_run = []
            while len(queue) > 0 and len(pairs_to_run) < self._num_cpus:
                pairs_to_run.append([queue.pop(0), 1])
            # Share out any remaining CPUs
            num_cpus_left = self._num_cpus - len(pairs_to_run)
            while num_cpus_left > 0:
                for i, _ in enumerate(pairs_to_run):
                    pairs_to_run[i][1] += 1
                    num_cpus_left -= 1
                    if num_cpus_left == 0: break
            
        return pairs_to_run
            

            
    def __hmm_command(self, input_pipe, pairs_to_run):
        r"""INTERNAL method for getting cmdline for running a batch of HMMs.
        
        Parameters
        ----------
        input_pipe: as hmmsearch
        pairs_to_run: list
            list with 2 members: (1) list of hmm and output file, (2) number of 
            CPUs to use when searching
            
        Returns
        -------
        A string command to be run with bash
        """
        element = pairs_to_run.pop()
        hmmsearch_cmd = self._individual_hmm_command(element[0][0],
                                                      element[0][1],
                                                      element[1])
        while len(pairs_to_run) > 0:
            element = pairs_to_run.pop()
            hmmsearch_cmd = "tee >(%s) | %s" % (self._individual_hmm_command(element[0][0],
                                                                              element[0][1],
                                                                              element[1]),
                                                hmmsearch_cmd)
            
        # Run the actual command
        hmmsearch_cmd = "%s | %s" % (input_pipe, hmmsearch_cmd)
        return hmmsearch_cmd
            
    def _individual_hmm_command(self, hmm, output_file, num_cpus):
        return "hmmsearch %s --cpu %s -o /dev/null --noali --domtblout %s %s -" % (self._extra_args,
                                                                         num_cpus,
                                                                         output_file,
                                                                         hmm)
        
class NhmmerSearcher(HmmSearcher):
    r"""Runs nhmmer given one or many HMMs in a scalable and fast way""" 
            
    def _individual_hmm_command(self, hmm, output_file, num_cpus):
        return "nhmmer %s --cpu %s -o /dev/null --noali --tblout %s %s -" % (self._extra_args,
                                                                   num_cpus,
                                                                   output_file,
                                                                   hmm)
        