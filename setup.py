from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('graftm/version.py').read()) # loads __version__

setup(name='graftm',
      version=__version__,
      author='Joel Boyd, Ben Woodcroft',
      description='GraftM is a pipeline used for identifying and classifying marker gene reads from metagenomic datasets',
      long_description=readme,
      license='GPL3+',
      keywords="",
      packages=find_packages(exclude='docs'),
      install_requires=('biopython >=1.64',
                        'seqmagick >=0.5.0',
                        'scikit-bio ==0.4.2',
                        'subprocess32 >=3.2.6',
                        'biom-format >=2.1.4',
                        'extern >=0.0.4',
                        'taxtastic >=0.5.4',
                        'tempdir >=0.6'),
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      url='http://geronimp.github.io/graftM',
      scripts=['bin/graftM'],
      data_files=[
          ('share', ['share/18S.hmm']),
      ],
)
