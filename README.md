# POSSEL a POSitive SELection pipeline


# Installation of Possel pipleine

## Instalation of requried Python libraries

### First create a conda environement possel
conda create -n possel
conda activate possel
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict

### install python 3.10 for conpatibility
conda install python==3.10

### installation of biopython (1.83 or higher)
pip install biopython==1.82 (avoid 1.83 as it has some issues with alignment)

### install python and basic datscience libraries
conda install pandas scipy statsmodels

### instalation of ete3 
pip install ete3

## instalation of the NextFLow workflow manager
conda install nextflow

## instalation of bioperl (it can get tricky sometime)
cpan install Bio
PERL5LIB=/Users/dthybert/perl5/lib/perl5
export PERL5LIB

## Instalation of extrnal tools
External tool:
 -MAFFT : 
 -zorro (includes FastTree)
 -pal2nal: https://www.bork.embl.de/pal2nal/#Download
 -raxml-ng: https://github.com/amkozlov/raxml-ng
 -hyphy : use the conda installation :conda install -c bioconda hyphiy
 -paml: https://github.com/abacus-gene/paml 
