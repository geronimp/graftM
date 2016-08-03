import json
import os
import shutil
import logging
import extern

from graftm.getaxnseq import Getaxnseq

class InsufficientGraftMPackageException(Exception): pass

class GraftMPackage:
    '''Package to represent a GraftM package. To start using a package, run

    pkg = GraftMPackage.acquire('/path/to/a.gpkg')

    and then retrieve values with e.g.

    pkg.alignment_hmm_path()
      #=> '/path/to/a.gpkg/the.hmm'
    '''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    # The key names are unlikely to change across package format versions,
    # so store them here in the superclass
    DIAMOND_DATABASE_KEY = 'diamond_database'
    VERSION_KEY = 'graftm_package_version'
    ALIGNMENT_HMM_KEY = 'align_hmm'
    SEARCH_HMM_KEY = "search_hmms"
    REFERENCE_PACKAGE_KEY = "refpkg"
    HMM_TRUSTED_CUTOFF_KEY = "trusted_cutoff"
    RANGE_KEY = "range"
    UNALIGNED_SEQUENCE_DATABASE_KEY = "unaligned_sequence_database"
    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    _CURRENT_VERSION = 3

    _REQUIRED_KEYS = {'2': [
                             VERSION_KEY,
                             ALIGNMENT_HMM_KEY,
                             SEARCH_HMM_KEY,
                             REFERENCE_PACKAGE_KEY,
                             HMM_TRUSTED_CUTOFF_KEY,
                             RANGE_KEY
                             ],
                      '3': [
                             VERSION_KEY,
                             ALIGNMENT_HMM_KEY,
                             SEARCH_HMM_KEY,
                             REFERENCE_PACKAGE_KEY,
                             HMM_TRUSTED_CUTOFF_KEY,
                             RANGE_KEY,
                             UNALIGNED_SEQUENCE_DATABASE_KEY
                             ]
                      }


    @staticmethod
    def acquire(graftm_package_path):
        '''Acquire a new graftm Package

        Parameters
        ----------
        graftm_output_path: str
            path to base directory of graftm
        '''
        
        contents_hash = json.load(
                                   open(
                                        os.path.join(
                                                     graftm_package_path,
                                                     GraftMPackage._CONTENTS_FILE_NAME
                                                     ),
                                         )
                                   )
        
        
        v=contents_hash[GraftMPackage.VERSION_KEY]
        logging.debug("Loading version %i GraftM package: %s" % (v, graftm_package_path))
        if v == 2:
            pkg = GraftMPackageVersion2()
        elif v == 3:
            pkg = GraftMPackageVersion3()
        else:
            raise InsufficientGraftMPackageException("Bad version: %s" % v)
        
        pkg._contents_hash = contents_hash
        pkg._base_directory = graftm_package_path
        # check we are at current version otherwise choke
        pkg.check_universal_keys(v)
        pkg.check_required_keys(GraftMPackage._REQUIRED_KEYS[str(v)])
        return pkg

    def check_universal_keys(self, version):
        h = self._contents_hash
        try:
            v = h[GraftMPackage.VERSION_KEY]
        except KeyError:
            raise InsufficientGraftMPackageException("No version information in graftm package")
        if v != version:
            raise InsufficientGraftMPackageException("Bad version: %s" % v)

    def check_required_keys(self, required_keys):
        '''raise InsufficientGraftMPackageException if this package does not
        conform to the standard of the given package'''
        h = self._contents_hash
        for key in required_keys:
            if key not in h:
                raise InsufficientGraftMPackageException("Package missing key %s" % key)

    def __getitem__(self, key):
        '''Return the value of the given key from the contents file'''
        return self.contents_hash[key]

    def contents_file_path(self):
        return os.path.join(self._base_directory, GraftMPackage._CONTENTS_FILE_NAME)

class GraftMPackageVersion2(GraftMPackage):
    version = 2

    def diamond_database_path(self):
        if self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY]:
            return os.path.join(self._base_directory,
                                self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY])
        else:
            return None
        
    def search_hmm_paths(self):
        return [os.path.join(self._base_directory, x) for x in
                self._contents_hash[GraftMPackage.SEARCH_HMM_KEY]]

    def alignment_hmm_path(self):
        return os.path.join(self._base_directory,
                            self._contents_hash[GraftMPackage.ALIGNMENT_HMM_KEY])
        
    def alignment_fasta_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['aln_fasta'])

    def reference_package_path(self):
        return os.path.join(self._base_directory,
                            self._contents_hash[GraftMPackage.REFERENCE_PACKAGE_KEY])

    def use_hmm_trusted_cutoff(self):
        return self._contents_hash[GraftMPackage.HMM_TRUSTED_CUTOFF_KEY]

    def maximum_range(self):
        return self._contents_hash[GraftMPackage.RANGE_KEY]

    def _refpkg_contents(self):
        return json.loads(open(os.path.join(self.reference_package_path(), 'CONTENTS.json')).read())

    def taxtastic_seqinfo_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['seq_info'])

    def taxtastic_taxonomy_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['taxonomy'])
        
    def reference_package_tree_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['tree'])

    @staticmethod
    def compile(output_package_path, refpkg_path, hmm_path, diamond_database_file, max_range, 
                trusted_cutoff=False, search_hmm_files=None):
        '''Create a new GraftM package with the given inputs. Any files
        specified as parameters are copied into the final package so can
        be removed after calling this function.

        Parameters
        ----------
        output_package_path: str
            path to the package being created (must not exist)
        refpkg_path: str
            path to pplacer reference package
        hmm_path: str
            path to the align HMM. Used as the search HMM if search_hmm_files
            is None
        diamond_database_file: str
            path to diamond DB file, or None for nucleotide packages
        max_rage: str
            as per maximum_range()
        trusted_cutoff: boolean
            set TC in search HMM

        search_hmm_files: list of str or None
            use these HMMs for search instead of the hmm_path. All
            basenames of these paths must be unique, and not the same as
            hmm_path.
            
        Returns
        -------
        Nothing
        '''

        if os.path.exists(output_package_path):
            raise Exception("Not writing new GraftM package to already existing file/directory with name %s" % output_package_path)
        os.mkdir(output_package_path)

        hmm_file_in_gpkg = os.path.basename(hmm_path)
        shutil.copyfile(hmm_path, os.path.join(output_package_path, hmm_file_in_gpkg))

        if diamond_database_file:
            diamond_database_file_in_gpkg = os.path.basename(diamond_database_file)
            shutil.copyfile(diamond_database_file, os.path.join(output_package_path, diamond_database_file_in_gpkg))

        refpkg_in_gpkg = os.path.basename(refpkg_path)
        shutil.copytree(refpkg_path, os.path.join(output_package_path, refpkg_in_gpkg))

        
        if search_hmm_files:
            search_hmm_files_in_gpkg_names =\
                [os.path.basename(path) for path in search_hmm_files]
            if hmm_file_in_gpkg in search_hmm_files_in_gpkg_names:
                raise Exception("Search and align HMMs must all have different basenames to create a new gpkg")
            if len(set(search_hmm_files_in_gpkg_names)) != len(search_hmm_files_in_gpkg_names):
                raise Exception("Search HMMs must have different basenames to create a new gpkg")
            for i, search_hmm in enumerate(search_hmm_files):
                shutil.copyfile(search_hmm, os.path.join(output_package_path, search_hmm_files_in_gpkg_names[i]))
        else:
            search_hmm_files_in_gpkg_names = [hmm_file_in_gpkg]
        
        contents = {GraftMPackage.VERSION_KEY: GraftMPackageVersion2.version,
                    GraftMPackage.ALIGNMENT_HMM_KEY: hmm_file_in_gpkg,
                    GraftMPackage.SEARCH_HMM_KEY: search_hmm_files_in_gpkg_names,
                    GraftMPackage.REFERENCE_PACKAGE_KEY: refpkg_in_gpkg,
                    GraftMPackage.HMM_TRUSTED_CUTOFF_KEY: trusted_cutoff,
                    GraftMPackage.RANGE_KEY: max_range}
        if diamond_database_file:
            contents[GraftMPackage.DIAMOND_DATABASE_KEY] = diamond_database_file_in_gpkg

        json.dump(contents, open(os.path.join(output_package_path, GraftMPackage._CONTENTS_FILE_NAME), 'w'))

    
class GraftMPackageVersion3(GraftMPackageVersion2):
    
    version = 3
    
    def unaligned_sequence_database_path(self):
        if self._contents_hash[GraftMPackage.UNALIGNED_SEQUENCE_DATABASE_KEY]:
            return os.path.join(self._base_directory,
                                self._contents_hash[GraftMPackage.UNALIGNED_SEQUENCE_DATABASE_KEY])
        else:
            return None

    def create_diamond_db(self):
        '''Create a diamond database from the unaligned sequences in this package.

        Returns
        -------
        path to the created diamond db e.g. 'my_sequences.dmnd'
        '''
        base = self.unaligned_sequence_database_path()
        cmd = "diamond makedb --in '%s' -d '%s'" % (self.unaligned_sequence_database_path(), base)
        extern.run(cmd)
        diamondb = '%s.dmnd' % base
        # Mostly this moves a file to it's current location because Create
        # follows this same logic, but there's a specially crafted
        # test/data/mcrA.gpkg which is slightly different.
        os.rename(diamondb, self.diamond_database_path())
        return diamondb

    def taxonomy_hash(self):
        '''Read in the taxonomy and return as a hash of name: taxonomy,
        where taxonomy is an array of strings.'''
        gtns = Getaxnseq()
        with open(self.taxtastic_taxonomy_path()) as tax:
            with open(self.taxtastic_seqinfo_path()) as seqinfo:
                return gtns.read_taxtastic_taxonomy_and_seqinfo(
                    tax, seqinfo)
        
    @staticmethod
    def graftm_package_is_protein(graftm_package):
        '''Return true if this package is an Amino Acid alignment package, otherwise
        False i.e. it is a nucleotide package. In general it is best to use
        'is_protein_package' instead.

        '''
        found = None
        with open(graftm_package.alignment_hmm_path()) as f:
            r = f.read().split("\n")
        for line in r:
            if line=='ALPH  DNA':
                found = False
                break
            elif line=='ALPH  amino':
                found = True
                break
        if found is None:
            raise Exception("Unable to determine whether the HMM was amino acid or dna")
        return found

    def is_protein_package(self):
        '''Return true if this package is an Amino Acid alignment package, otherwise
        False i.e. it is a nucleotide package. Cache the result for speed.

        '''
        if not hasattr(self, '_is_protein_package'):
            self._is_protein_package = GraftMPackageVersion3.graftm_package_is_protein(self)
        return self._is_protein_package

        
    @staticmethod
    def compile(output_package_path, refpkg_path, hmm_path, diamond_database_file, 
                max_range, unaligned_sequence_database=None, 
                trusted_cutoff=False, search_hmm_files=None):
        '''Create a new GraftM package with the given inputs. Any files
        specified as parameters are copied into the final package so can
        be removed after calling this function.
        
        Parameters
        ----------
        output_package_path: str
            path to the package being created (must not exist)
        refpkg_path: str
            path to pplacer reference package
        hmm_path: str
            path to the search and align HMM
        diamond_database_file: str
            path to diamond DB file, or None for nucleotide packages
        max_rage: str
            as per maximum_range()
        unaligned_sequence_database: str
            path to unaligned sequence database
        trusted_cutoff: boolean
            set TC in search HMM
        search_hmm_files: list of str or None
            use these HMMs for search instead of the hmm_path. All
            basenames of these paths must be unique, and not the same as
            hmm_path.
            
        Returns
        -------
        Nothing
        '''
        
        if os.path.exists(output_package_path): 
            raise Exception("Not writing new GraftM package to already existing file/directory with name %s" % output_package_path)
        os.mkdir(output_package_path)
        
        hmm_file_in_gpkg = os.path.basename(hmm_path)
        shutil.copyfile(hmm_path, os.path.join(output_package_path, hmm_file_in_gpkg))
        # Copy unaligned sequence database into graftm package
        if unaligned_sequence_database:
            unaligned_sequence_database_in_gpkg = os.path.join(output_package_path, os.path.basename(unaligned_sequence_database))
            shutil.copyfile(unaligned_sequence_database, unaligned_sequence_database_in_gpkg)
        
        if diamond_database_file:
            diamond_database_file_in_gpkg = os.path.basename(diamond_database_file)
            shutil.copyfile(diamond_database_file, os.path.join(output_package_path, diamond_database_file_in_gpkg))
        else:
            diamond_database_file_in_gpkg = diamond_database_file
        refpkg_in_gpkg = os.path.basename(refpkg_path)
        shutil.copytree(refpkg_path, os.path.join(output_package_path, refpkg_in_gpkg))
        
        if search_hmm_files:
            search_hmm_files_in_gpkg_names =\
                [os.path.basename(path) for path in search_hmm_files]
            if hmm_file_in_gpkg in search_hmm_files_in_gpkg_names:
                raise Exception("Search and align HMMs must all have different basenames to create a new gpkg")
            if len(set(search_hmm_files_in_gpkg_names)) != len(search_hmm_files_in_gpkg_names):
                raise Exception("Search HMMs must have different basenames to create a new gpkg")
            for i, search_hmm in enumerate(search_hmm_files):
                shutil.copyfile(search_hmm, os.path.join(output_package_path, search_hmm_files_in_gpkg_names[i]))
        else:
            search_hmm_files_in_gpkg_names = [hmm_file_in_gpkg]
        
        contents = {GraftMPackage.VERSION_KEY: GraftMPackageVersion3.version,
                    GraftMPackage.ALIGNMENT_HMM_KEY: hmm_file_in_gpkg,
                    GraftMPackage.SEARCH_HMM_KEY: search_hmm_files_in_gpkg_names,
                    GraftMPackage.REFERENCE_PACKAGE_KEY: refpkg_in_gpkg,
                    GraftMPackage.HMM_TRUSTED_CUTOFF_KEY: trusted_cutoff,
                    GraftMPackage.RANGE_KEY: max_range,
                    GraftMPackage.UNALIGNED_SEQUENCE_DATABASE_KEY: (os.path.basename(unaligned_sequence_database)
                                                                    if unaligned_sequence_database else None),
                    GraftMPackage.DIAMOND_DATABASE_KEY: diamond_database_file_in_gpkg}
        
        json.dump(contents, open(os.path.join(output_package_path, GraftMPackage._CONTENTS_FILE_NAME), 'w'))
