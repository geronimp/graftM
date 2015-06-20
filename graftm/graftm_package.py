import json
import os

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
    VERSION_KEY = 'version'
    ALIGNMENT_HMM_KEY = 'aln_hmm'
    SEARCH_HMM_KEY = "search_hmm"
    REFERENCE_PACKAGE_KEY = "rfpkg"
    HMM_TRUSTED_CUTOFF_KEY = "TC"
    
    @staticmethod
    def acquire(graftm_package_path):
        '''Acquire a new graftm Package
        
        Parameters
        ----------
        graftm_output_path: str
            path to base directory of graftm
        '''
        pkg = GraftMPackageVersion2()
        
        pkg._base_directory = graftm_package_path
        pkg._contents_hash = json.load(open(os.path.join(graftm_package_path,
                                                        GraftMPackage._CONTENTS_FILE_NAME), 
                                             ))
        
        # check we are at current version otherwise choke
        pkg.check_universal_keys(2)
        pkg.check_required_keys(GraftMPackageVersion2._REQUIRED_KEYS)
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
                raise InsufficientGraftMPackageException("package missing key %s" % key)
            
    def __getitem__(self, key):
        '''Return the value of the given key from the contents file'''
        return self._contents_hash[key]
    
class GraftMPackageVersion2(GraftMPackage):
    version = 2
    
    _REQUIRED_KEYS = [
                     GraftMPackage.DIAMOND_DATABASE_KEY,
                     GraftMPackage.VERSION_KEY,
                     GraftMPackage.ALIGNMENT_HMM_KEY,
                     GraftMPackage.SEARCH_HMM_KEY,
                     GraftMPackage.REFERENCE_PACKAGE_KEY,
                     GraftMPackage.HMM_TRUSTED_CUTOFF_KEY
                     ]
    
    def diamond_database_path(self):
        return os.path.join(self._base_directory, 
                            self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY])
    
    def search_hmm_paths(self):
        return [os.path.join(self._base_directory, x) for x in
                self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY]]
        
    def alignment_hmm_path(self):
        return os.path.join(self._base_directory, 
                            self._contents_hash[GraftMPackage.ALIGNMENT_HMM_KEY])
        
    def reference_package_path(self):
        return os.path.join(self._base_directory, 
                            self._contents_hash[GraftMPackage.REFERENCE_PACKAGE_KEY])
        
    def use_hmm_trusted_cutoff(self):
        return self._contents_hash[GraftMPackage.HMM_TRUSTED_CUTOFF_KEY]
        
        
    
    

            