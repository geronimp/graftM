class MalformedGreenGenesTaxonomyException(Exception):  pass
class DuplicateTaxonomyException(Exception):  pass

class GreenGenesTaxonomy:
    # A hash of name to array of taxonomy strings
    taxonomy = {}

    def __init__(self, taxonomy={}):
        self.taxonomy = taxonomy

    @staticmethod
    def read(input_taxonomy_io):
        '''Parse in a taxonomy from the given I/O.

        Parameters
        ----------
        input_taxonomy_io: io
            an open io object that is iterable. This stream is neither opened
            nor closed by this method

        Returns
        -------
        A GreenGenesTaxonomy instance with the taxonomy parsed in

        Raises
        ------
        DuplicateTaxonomyException: if there is duplicate taxonomy entries
        MalformedGreenGenesTaxonomyException: if there is something else amiss'''
        tax = {}
        for line in input_taxonomy_io:
            if len(line.strip()) == 0: continue #ignore empty lines
            splits = line.split("\t")
            if len(splits) != 2:
                raise MalformedGreenGenesTaxonomyException("Unexpected number of tab-separated fields found in taxonomy file, on line %s" % line)
            name = splits[0].strip()
            taxonomy = [t.strip() for t in splits[1].split(';')]
            while len(taxonomy) > 0 and taxonomy[-1] == '':
                taxonomy = taxonomy[:-1]
            for lineage in taxonomy:
                if lineage == '':
                    raise MalformedGreenGenesTaxonomyException("Encountered a taxonomy string with the middle of the taxonomy string missing: %s" % line)

            if name in tax:
                raise DuplicateTaxonomyException("Duplicate definition of taxonomy for %s" % name)
            tax[name] = taxonomy
        return GreenGenesTaxonomy(tax)

    @staticmethod
    def read_file(input_filename):
        '''Like read() except uses a file rather than an IO stream, for convenience'''
        with open(input_filename) as f:
            g = GreenGenesTaxonomy.read(f)
        return g

    def write(self, output_io):
        '''Write a taxonomy to an open stream out in GG format. Code calling this
        function must open and close the io object.'''
        for name, tax in self.taxonomy.items():
            output_io.write("%s\t%s\n" % (name, '; '.join(tax)))
