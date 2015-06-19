import json
import os

class InsufficientGraftMPaackageException(Exception): pass

class GraftMPackage:
    '''Package to represent a graftm package. To start using a package, run
    pkg = GraftMPackage.acquire, and then retrieve values with e.g. 
    pkg[GraftMPackage.VERSION_KEY]'''
     
    DIAMOND_DATABASE_KEY = 'diamond_database'
    VERSION_KEY = 'version'
    ALIGNMENT_HMM_KEY = 'aln_hmm'
    SEARCH_HMM_KEY = "search_hmm"
    REFERENCE_PACKAGE_KEY = "rfpkg"
    HMM_TRUSTED_CUTOFF_KEY = "TC"
    TAX_INFO_KEY = "tax_info"
    _CONTENTS_FILE_NAME = 'CONTENTS.json'
    
    _CURRENT_VERSION = 2
    
    _REQUIRED_KEYS = {'2': [
                     DIAMOND_DATABASE_KEY,
                     VERSION_KEY,
                     ALIGNMENT_HMM_KEY,
                     SEARCH_HMM_KEY,
                     REFERENCE_PACKAGE_KEY,
                     HMM_TRUSTED_CUTOFF_KEY,
                     TAX_INFO_KEY
                     ]}
    
    @staticmethod
    def acquire(graftm_package_path):
        '''Acquire a new graftm Package
        
        Parameters
        ----------
        graftm_output_path: str
            path to base directory of graftm
        '''
        pkg = GraftMPackage()
        
        
        pkg._base_directory = graftm_package_path
        pkg._contents_hash = json.load(open(os.path.join(graftm_package_path,
                                                        GraftMPackage._CONTENTS_FILE_NAME), 
                                             ))
        
        # check we are at current version otherwise choke
        
        pkg.check_required_keys(pkg._CURRENT_VERSION)
        return pkg
        
        
    def check_required_keys(self, version):
        '''raise InsufficientGraftMPaackageException if this package does not
        conform to the standard of the given package'''
        #check version key
        h = self._contents_hash
        v=(h[self.VERSION_KEY] if self.VERSION_KEY in h else None)
        if not v:
            raise InsufficientGraftMPaackageException("No version information in graftm package")
        if v != version:
            raise InsufficientGraftMPaackageException("Bad version: %s" % v)
        
        for key in self._REQUIRED_KEYS[str(version)]:
            if key not in h:
                raise InsufficientGraftMPaackageException("package missing key %s" % key)
            
    def __getitem__(self, key):
        '''Return the value of the given key from the contents file'''
        return self._contents_hash[key]
            