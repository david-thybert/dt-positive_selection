import typing
import sys
import os
import argparse 
from Bio import SeqIO



def main(pep_align: str, nucleotide: str, outBase: str)-> None:
    pass


parser = argparse.ArgumentParser(description='The script check that the nucleotide unaligned file correspsond to the aligned peptide file')
parser.add_argument('--pepal',type=str, help='peptide alignemnt file')
parser.add_argument('--nuc',type=str, help='nucleotide unalign filee')
parser.add_argument('--out', type=str, help='path to the outdirectory')
args = parser.parse_args()

main(args.pepal, args.nuc, args.out)