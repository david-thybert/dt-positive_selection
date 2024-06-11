import argparse 
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
from Bio import AlignIO


def main(codeml_command:str, ctl_file:str)->None:

    command = f"{codeml_command} {ctl_file}"
    out_process = subprocess.run(command, shell=True, capture_output=True, text=True)
    if out_process.stderr.split("\n")[1] != 'Error: end of tree file..':
        print(f"error running codeml with {ctl_file}")
        raise Exception("error running codeml", ctl_file)


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--codeml',type=str, help='multiple alignment in fasta format')
parser.add_argument('--ctl', type=str, help='threshold used to remove unconfident alignment')


args = parser.parse_args()

main(args.codeml, args.ctl)
