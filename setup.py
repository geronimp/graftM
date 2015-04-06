from disutils.core import setup

setup(
    name='GraftM',
    version='0.4.2',
    author='Joel Boyd, Ben Woodcroft',
    author_email='joel.boyd@uqconnect.edu.au',
    packages=['graftm', 'bin'],
    url='http://pypi.python.org/pypi/graftm/',
    license='GPL3',
    description='derives community composition from short read data',
    long_description=open('README.txt').read(),
    install_requires=[])
