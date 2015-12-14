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
* diamond v. >= 0.7.9 (https://github.com/bbuchfink/diamond)

To create new GraftM packages, you'll also need
* FastTreeMP (http://www.microbesonline.org/fasttree/)

#### Docker images
Versions of graftM on pip now have matching docker images as of GraftM v0.9.2. GraftM docker images are portable containers that contain the graftM code and all its python and non-python dependancies, allowing GraftM to be run on any platform with docker installed. Details on how to download and run a GraftM docker image can be found on the [graftm-docker](https://github.com/geronimp/graftM-docker) GitHub page or the [docker hub page](https://hub.docker.com/u/geronimp/).

### Manual
A [manual](https://github.com/geronimp/graftM/wiki) is available in the form of the wiki here on GitHub.

### GraftM packages
We have a starter pack of graftM packages available including:

* 16S rRNA packages
* 15 single copy ribosomal protein marker genes
* The methanogenesis marker mcrA

All are available at the [GraftM package database store](https://drive.google.com/open?id=0BwJ4AwdqUiTzfndmRXowX3MydkM5bG1PYmxRUjNmMUNnazdFaUJaWjJFSkh1UEFDSkpReU0).

Once you have downloaded the package you want, just decompress it as follows:

```
tar -xvzf my.tar.gz
```
And you should be good to go!



### Example
As an example, we'll use GraftM to classify a single 16S sequence from GreenGenes. Saving the example file as `/tmp/eg.fasta` with the following contents:
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

### Contact
If you have any further comments, complaints or recomendations about GraftM, drop an email to the [SupportM](https://groups.google.com/forum/?hl=en#!forum/supportm) public help forum.
Software by [Joel A. Boyd](http://ecogenomic.org/users/joel-boyd) (@geronimp) and [Ben J. Woodcroft](http://www.ecogenomic.org/users/ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics](http://ecogenomic.org)
Released under GPL3 - see LICENCE.txt for licencing details
