import argparse

from Bio import AlignIO
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter



"""
This script convert a multiple sequence alignemnt in fasta to phylip format

command:

    python fasta2phy.py --mult_fasta <fasta file> --out_phy <alignemnt file .phy>

"""

def main(mult_fasta:str, mult_phy:str)-> None:
    """
    This script convert a multiple alignment in fsta format 
    to phylip format

    :param mult_fasta: path to the mulitple alignment filw in fasta format
    :param mult_phy: path to the outfile in phylipe format 
    """
    with open(mult_fasta) as input_handle:
        alignment = AlignIO.read(input_handle, 'fasta')
    with open(mult_phy, 'w') as output_handle:  
        for seq in alignment:
            print(seq.id)
            seq.id = seq.id + " "
        SequentialPhylipWriter(output_handle).write_alignment(alignment, id_width=50)



########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script converting a multiple sequence alignment from fast to phylip format')

parser.add_argument('--mult_fasta',type=str, help='multiple alignment in fasta format')
parser.add_argument('--out_phy', type=str, help='name of the phyle in phylip format')

args = parser.parse_args()
main(args.mult_fasta, args.out_phy)