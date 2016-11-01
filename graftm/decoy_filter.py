import logging
from graftm.diamond import Diamond
from graftm.unpack_sequences import UnpackRawReads
from graftm.sequence_search_results import SequenceSearchResult
from graftm.sequence_extractor import SequenceExtractor

class DecoyFilter:
    def __init__(self, proper_hits_diamond, decoy_diamond=None):
        '''Create a new DecoyFilter.

        Parameters
        ----------
        proper_hits_diamond: Diamond
            Diamond object containing real sequences.
        decoy_diamond: Diamond or None.
            Diamond object with db containing decoy sequences. If None, no
            searching against a decoy database is carried out.

        '''
        self._decoy_diamond = decoy_diamond
        self._proper_hits_diamond = proper_hits_diamond
        
    def filter(self, candidate_sequences_fasta_path, filtered_output_fasta_path):
        '''Filter the fasta file by only keeping sequences that hit a proper sequence
        better than the decoy database.
            
        Parameters
        ----------
        candidate_sequences_fasta_path: str
            path to candidate sequences fasta file i.e. sequences to be filtered.
        filtered_output_fasta_path: str
            path to file which, after this method is run, will contain only
            those sequences which pass the filter.
            
        Returns
        -------
        False if no sequences remain after filtering, else True.
        '''
        # Run query sequences against the proper database
        seq_ids_and_bitscores = {}
        logging.debug("Running diamond against the non-decoy sequences")
        pd = self._proper_hits_diamond.run(
            candidate_sequences_fasta_path,
            UnpackRawReads.PROTEIN_SEQUENCE_TYPE)
        for res in pd.each([SequenceSearchResult.QUERY_ID_FIELD,
                            SequenceSearchResult.ALIGNMENT_BIT_SCORE]):
            seq = res[0]
            score = float(res[1])
            # Possible a single sequence gets 2 split up hits (maybe), so take
            # the highest bitscore.
            if seq not in seq_ids_and_bitscores or seq_ids_and_bitscores[seq] < score:
                seq_ids_and_bitscores[seq] = score

        num_before_decoy_removal = len(seq_ids_and_bitscores)
        logging.info("Found %i sequences which hit the non-decoy sequences" %\
                     num_before_decoy_removal)

        if self._decoy_diamond is None:
            logging.debug("Not running against the decoy database")
        else:
            # Run the query sequences against the decoy database, removing from the
            # list any sequences which hit better the decoy DB.
            logging.debug("Running diamond against decoy sequences")
            pd = self._decoy_diamond.run(
                candidate_sequences_fasta_path,
                UnpackRawReads.PROTEIN_SEQUENCE_TYPE)
            for res in pd.each([SequenceSearchResult.QUERY_ID_FIELD,
                                SequenceSearchResult.ALIGNMENT_BIT_SCORE]):
                seq = res[0]
                score = float(res[1])
                if seq in seq_ids_and_bitscores and seq_ids_and_bitscores[seq] < score:
                    logging.debug("Removing sequence with better hit to the decoy database: %s" % seq)
                    del seq_ids_and_bitscores[seq]

            logging.info("Removed %i"
                         " sequences which hit the decoy sequences better"
                         " than the non-decoy sequences" %\
                         (num_before_decoy_removal-len(seq_ids_and_bitscores)))
            if len(seq_ids_and_bitscores) == 0:
                return False

        # Extract the found sequences into the output file
        logging.debug("Extracting query sequences")
        SequenceExtractor().extract(seq_ids_and_bitscores.keys(),
                                    candidate_sequences_fasta_path,
                                    filtered_output_fasta_path)
        return True
