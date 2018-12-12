# GraftM
GraftM is a tool for finding genes of interest in metagenomes, metatranscriptomes, and whole genomes.

Using modular gene packages, GraftM will search the provided sequences using hmmsearch (HMMER) and place the identified sequences into a pre-constructed phylogenetic tree. The provides fast, phylogenetically informed community profiles and genome annotations. GraftM provides tools to:
* Create and update custom gene packages to use with GraftM
* Decorate trees, and of course..
* Analyse sequence datasets using these GraftM packages
Head over to the [GraftM page](http://geronimp.github.io/graftM/) for more general information.

### Installation
#### pip installation
GraftM can be installed via pip:
```
pip install graftm
```
However, to use all features of GraftM a few extra binary applications are required:
* orfm v. >= 0.2.0 (https://github.com/wwood/OrfM)
* hmmer v. >= 3.1b1 (http://hmmer.janelia.org/)
* fxtract (https://github.com/ctSkennerton/fxtract)
* pplacer v. >= 2.6.32 (http://matsen.fhcrc.org/pplacer/)
* krona v. >= 2.4 (http://sourceforge.net/p/krona/home/krona/)
* mafft v. >= 7.22 (http://mafft.cbrc.jp/)
* diamond v. >= 0.9 (https://github.com/bbuchfink/diamond) The version of diamond used must be matched with that used to generate gpkgs, and to run the tests a specific version is required.

To create new GraftM packages, you'll also need
* FastTreeMP (http://www.microbesonline.org/fasttree/)

#### Docker images
Versions of graftM on pip now have matching docker images as of GraftM v0.9.2. GraftM docker images are portable containers that contain the graftM code and all its python and non-python dependancies, allowing GraftM to be run on any platform with docker installed. Details on how to download and run a GraftM docker image can be found on the [graftm-docker](https://github.com/geronimp/graftM-docker) GitHub page or the [docker hub page](https://hub.docker.com/u/geronimp/).

#### GNU Guix
GraftM and all its Python and non-Python dependencies can be installed using the ACE package repository. After installing [GNU Guix](https://www.gnu.org/software/guix/)
```
git clone https://github.com/Ecogenomics/ace-guix
GUIX_PACKAGE_PATH=ace-guix guix package --install graftm
```

#### Conda
Although not officially supported as an installation method, GraftM can be installed through [conda / bioconda](https://anaconda.org/bioconda/graftm).

### Manual
A [manual](https://github.com/geronimp/graftM/wiki) is available in the form of the wiki here on GitHub.

### GraftM packages
We have a starter pack of graftM packages available including:

* 16S rRNA packages
* 15 single copy ribosomal protein marker genes
* The methanogenesis marker mcrA

All are available at the [ACE ftp GraftM package database store](https://data.ace.uq.edu.au/public/graftm).

Once you have downloaded the package you want, just decompress it as follows:

```
tar -xvzf my.tar.gz
```
And you should be good to go!



### Example: 'graft'
We'll use GraftM to classify a single 16S sequence from GreenGenes. Save the example file as `/tmp/eg.fasta` with the following contents:
```
>229854
GAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGCAGCATGACTTAGCTTGCT
AAGTTGATGGCGAGTGGCGAACGGGTGAGTAACGCGTAGGAATATGCCTTAAAGAGGGGGACAACTTGGGGAAACTCAAG
CTAATACCGCATAAACTCTTCGGAGAAAAGCTGGGGACTTTCGAGCCTGGCGCTTTAAGATTAGCCTGCGTCCGATTAGC
```
Then we can use GraftM's 61% OTU clustered GraftM package to detect and classify this sequence. Running graftM might look something like this:
```
$ graftM graft --forward /tmp/eg.fasta --graftm_package /path/to/4.01.2013_08_greengenes_61_otus.gpkg/ --output_directory eg.graftm

                             GraftM 0.9.2

                                GRAFT

                       Joel Boyd, Ben Woodcroft

                                                         __/__
                                                  ______|
          _- - _                         ________|      |_____/
           - -            -             |        |____/_
           - _     >>>>  -   >>>>   ____|
          - _-  -         -             |      ______
             - _                        |_____|
           -                                  |______

12/02/2015 09:52:06 AM INFO: Working on eg
12/02/2015 09:52:06 AM INFO: 1 read(s) detected
12/02/2015 09:52:06 AM INFO: aligning reads to reference package database
12/02/2015 09:52:06 AM INFO: Filtered 0 short sequences from the alignment
12/02/2015 09:52:06 AM INFO: 1 sequences remaining
12/02/2015 09:52:06 AM INFO: Placing reads into phylogenetic tree
12/02/2015 09:52:07 AM INFO: Placements finished
12/02/2015 09:52:07 AM INFO: Reading classifications
12/02/2015 09:52:07 AM INFO: Reads classified
12/02/2015 09:52:07 AM INFO: Writing summary table
12/02/2015 09:52:07 AM INFO: Writing biom file
12/02/2015 09:52:07 AM INFO: Building summary krona plot
12/02/2015 09:52:07 AM INFO: Cleaning up
12/02/2015 09:52:07 AM INFO: Done, thanks for using graftM!
```
This creates a folder `eg.graftm` which contains the results.

### Example: 'create'

We encourage you to create your own GraftM packages that suit your own analysis
needs (e.g. new genes, more representatives, improved annotation of phylogenetic
trees, etc.). Here we will create a 16S package out from a highly dereplicated
collection provided by GreenGenes. The data here are provided in the
[example_data/create](https://github.com/geronimp/graftM/tree/master/example_data/create)
folder of GraftM.
```
$ graftM create --output /tmp/my.gpkg --sequences 61_otus.fasta --taxonomy 61_otu_taxonomy.txt
                         
                             GraftM 0.9.5

                            CREATE

                   Joel Boyd, Ben Woodcroft

                                                    /
              >a                                   /
              -------------                       /
              >b                        |        |
              --------          >>>     |  GPKG  |
              >c                        |________|
              ----------

08/02/2016 09:18:48 AM WARNING: Deleting previous directory /tmp/my.gpkg
08/02/2016 09:18:48 AM INFO: Building gpkg for /tmp/my.gpkg
08/02/2016 09:18:48 AM INFO: Building seqinfo and taxonomy file from input taxonomy
08/02/2016 09:18:48 AM INFO: Checking for duplicate sequences
08/02/2016 09:18:48 AM INFO: Aligning sequences to create aligned FASTA file
08/02/2016 09:19:08 AM INFO: Building HMM from alignment
08/02/2016 09:19:09 AM INFO: Filtered 0 short sequences from the alignment
08/02/2016 09:19:09 AM INFO: 22 sequences remaining
08/02/2016 09:19:09 AM INFO: Checking for incorrect or fragmented reads
08/02/2016 09:19:09 AM INFO: Removing 0 sequences from the search HMM that are redundant at the 6 rank in the taxonomy file
08/02/2016 09:19:26 AM INFO: Building HMM from alignment
08/02/2016 09:19:28 AM INFO: Filtered 0 short sequences from the alignment
08/02/2016 09:19:28 AM INFO: 22 sequences remaining
08/02/2016 09:19:28 AM WARNING: Found a non-standard character in the sequence of 4363260: e.g. 'W'
08/02/2016 09:19:28 AM INFO: Deduplicating sequences
08/02/2016 09:19:28 AM INFO: Removed 0 sequences as duplicates, leaving 22 non-identical sequences
08/02/2016 09:19:28 AM INFO: Building tree
08/02/2016 09:19:29 AM INFO: Building seqinfo and taxonomy file from input taxonomy
08/02/2016 09:19:29 AM INFO: Creating reference package
08/02/2016 09:19:29 AM INFO: Attempting to run taxit create with rerooting capabilities
08/02/2016 09:19:29 AM INFO: Creating diamond database
08/02/2016 09:19:29 AM INFO: Compiling gpkg
08/02/2016 09:19:29 AM INFO: Cleaning up
08/02/2016 09:19:29 AM INFO: Testing gpkg package works
08/02/2016 09:19:37 AM INFO: Finished
```

Then, the output package `/tmp/my.gpkg` can be provided to 'graft'. There are
many optional arguments to 'create' which can be used to modify the process of
making gpkgs.

Happy grafting!

### Contact
If you have any further comments, complaints or recomendations about GraftM, drop an email to the [SupportM](https://groups.google.com/forum/?hl=en#!forum/supportm) public help forum.
Software by [Joel A. Boyd](http://ecogenomic.org/users/joel-boyd) (@geronimp) and [Ben J. Woodcroft](http://www.ecogenomic.org/users/ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics](http://ecogenomic.org)
Released under GPL3 - see LICENCE.txt for licencing details

### Citation
GraftM has been published - please cite us at

> **GraftM: a tool for scalable, phylogenetically informed classification of genes within metagenomes**.
> Joel A Boyd Ben J Woodcroft Gene W Tyson
> Nucleic Acids Research, Volume 46, Issue 10, 1 June 2018, Pages e59, https://doi.org/10.1093/nar/gky174

### License
GraftM is licensed under the GNU GPL v3+. See LICENSE.txt for further details. GraftM makes use of 18S rRNA sequences sourced from the [SILVA database](https://www.arb-silva.de), which employs a dual licensing model. See SILVA.LICENCE.txt for further details.
