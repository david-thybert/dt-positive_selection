import argparse 
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
from Bio import AlignIO

"""
This script launch codeml and check that there is not errors

command:

    python run_codeml.py --codeml <codeml command> --ctl <ctl config file>

"""

def main(codeml_command:str, ctl_file:str)->None:
    """
    The main funciton fo the script

    :param codeml_command: codeml commamd
    :param ctl_file: path to the configuration ctl file
    """
    command = f"{codeml_command} {ctl_file}"
    out_process = subprocess.run(command, shell=True, capture_output=True, text=True)
    if out_process.stderr.split("\n")[1] != 'Error: end of tree file..':
        print(f"error running codeml with {ctl_file}")
        raise Exception("error running codeml", ctl_file)


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='This cript run codeml')
parser.add_argument('--codeml',type=str, help='multiple alignment in fasta format')
parser.add_argument('--ctl', type=str, help='threshold used to remove unconfident alignment')

args = parser.parse_args()
main(args.codeml, args.ctl)
