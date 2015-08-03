import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('graftm/version.py').read()) # loads __version__

setup(name='graftm',
      version=__version__,
      author='Joel Boyd, Ben Woodcroft',
      description='GraftM is a pipeline used for identifying and classifying marker gene reads from metagenomic datasets',
      long_description=readme,
      license='see LICENSE.txt',
      keywords="",
      packages=find_packages(exclude='docs'),
      install_requires=('biopython ==1.64',
                        'seqmagick ==0.5.0',
                        'scikit-bio ==0.2.2',
                        'subprocess32 ==3.2.6',
                        'biom-format ==2.1.4',
                        'extern ==0.0.4'),
      url='http://geronimp.github.io/graftM',
      scripts=['bin/graftM']
)
