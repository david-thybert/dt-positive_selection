import argparse
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
from Bio import AlignIO
import pandas as pd

"""
This script identify with zorro region of a multiple sequence alignment that is not confident and remove them.

command:

    python zorro_wrapper.py --zorro <path to zorro> --tree_cmd <path to Fastree> --mult <fmultiple alignment fasta> --threshold <0. - 10.> --out <filtered alignment> 

"""


def run_zorro(mult_fasta: str, zorro_cmd: str, tree_cmd: str) -> list:
    """
    Wrapper around zorro and return the score for the alignment in each position

    :param mult_fasta: multiple alignment in fasta
    :param zorro_cmd: path  to zorro executable
    :param tree_cmd: tree tool needed by zorro
    :return: the list of score corresponding at each column of the alignment
    """
    result = []
    command = f"{zorro_cmd} -treeprog {tree_cmd} {mult_fasta}"
    out_process = subprocess.run(command, shell=True, capture_output=True)
    print(command)
    print(out_process.returncode)
    # if out_process.returncode != 1: # strangly zorro return 1 when process succesful
    #    raise Exception(f"issue running zorro with file {mult_fasta}")

    try:  # this try blcok is to compensate the non follwoing the norm of zorro with eror status
        lst_score = out_process.stdout.split()
        for score in lst_score:
            result.append(float(score))
    except Exception as e:
        print(f"issue running zorro with file {mult_fasta}")
        raise e
    return result


def load_align(mult_fasta: str) -> object:
    """
    Load multiple peptide sequence alignment

    :param mult_fasta:  multiple alignement file in fasta format
    :return: Alignment object
    """
    alignment = None
    with open(mult_fasta) as file_handler:
        alignment = AlignIO.read(file_handler, "fasta")
    return alignment


def get_confident_regions(scored_pos: list, threshold: float) -> list:
    """
    Get a list fo region fo confident alignment from zorro output.
    This function works at the peptide level

    :param confidence_pos: list fo confidence score per position ordered
    :param thereshold: threshold define a confident/non confident alignemnt
    """
    i = 0
    state = 0  # state 0 : belongs to a confident aligned position
    # state 1 : belong to a non confident aligned position
    start = 0
    end = -1
    confident = []
    while i < len(scored_pos):
        if scored_pos[i] < threshold:
            if state == 0:
                end = i
                if end != start:
                    confident.append([start, end])
                state = 1
            elif state == 1:
                pass
        else:
            if state == 0:
                pass
            elif state == 1:
                start = i
                state = 0
        i = i + 1

    if state == 0 and start < i and end == -1:
        confident.append([start, i])
    elif state == 0 and start < i and end <= i:
        confident.append([start, i])
    elif state == 0 and start < i and end >= i:
        confident.append([start, i])

    return confident


def filter_align_pep(alignment: object, confidence_pos: list) -> object:
    """
    The function keep the high confident column of the alignment in a
    new alignemnt object

    :param alignment: alignment object
    :param confidence_pos: list of confident regions in the fomrat of [[start:end],...]
    :return: concatenation of the confident regions in the order in an alignemnt aobject
    """
    filtered_alignement = None
    print(confidence_pos)
    if confidence_pos == []:
        return filtered_alignement
    for interval in confidence_pos:
        start = interval[0]
        end = interval[-1]
        if filtered_alignement is None:
            filtered_alignement = alignment[:, start:end]
        else:
            filtered_alignement = filtered_alignement + alignment[:, start:end]
    print(format(filtered_alignement, "fasta"))
    return filtered_alignement


def filter_align_nuc(alignment: object, confidence_pos: list) -> object:
    """
    The function keep the high confident column of the alignment in a
    new alignemnt object

    :param alignment: alignment object
    :param confidence_pos: list of confident regions in the fomrat of [[start:end],...] in aminoacid coordinate
    :return: concatenation of the confident regions in the order in an alignemnt aobject
    """
    filtered_alignement = None
    if confidence_pos == []:
        return filtered_alignement
    for interval in confidence_pos:
        start = interval[0] * 3
        end = interval[-1] * 3
        if filtered_alignement is None:
            filtered_alignement = alignment[:, start:end]
        else:
            filtered_alignement = filtered_alignement + alignment[:, start:end]
    print(format(filtered_alignement, "fasta"))
    return filtered_alignement

def load_score(guid_score:str)-> list:
    result = []
    with open(guid_score) as file_hand:
        for line in file_hand:
            if "#" in line:
                continue
            tab = line.split("\t")
            result.append(float(tab[1]))
    return result


def main(mult_pep_fasta: str, mult_nuc_fasta: str, guid_score: str, threshold: float, out_file: str) -> None:
    """
    This is the main function of the script

    :param mult_pep_fasta: multiple alignment in fasta format
    :param mult_nuc_fasta: path to zorro executable
    :param guid_score: command for tree tool
    :param threshold: threshold for filtering no confident alignment
    :param out_file: out_File storing the filtered alignment
    :return: nothing
    """
    # run zorro to score alignment column confidence
    scored_pos = load_score(guid_score)

    # identify highly confident position based on threshold
    confidence_pos = get_confident_regions(scored_pos, threshold)
    print(confidence_pos)

   # keep only highly confident regions from alignment
    print(mult_pep_fasta)
    alignment_pep = load_align(mult_pep_fasta)
    alignment_pep_filtered = filter_align_pep(alignment_pep, confidence_pos)

    # keep only highly confident regions from alignment
    alignment_nuc = load_align(mult_nuc_fasta)
    alignment_nuc_filtered = filter_align_nuc(alignment_nuc, confidence_pos)

    # save the filtered alignment
    with open(f"{out_file}.pep.filt", "w") as file_hanlder:
        file_hanlder.write(format(alignment_pep_filtered, "fasta"))

        # save the nuc filtered alignment
    with open(f"{out_file}.nuc.filt", "w") as file_hanlder:
        file_hanlder.write(format(alignment_nuc_filtered, "fasta"))

    # save the regions files
    regions_file = f"{out_file}.confident_reg"
    with open(regions_file, "w") as file_handler:
        print(confidence_pos)
        for pos in confidence_pos:
            file_handler.write(f"{pos[0]}\t{pos[-1]}\n")


########################################################################################
########### Main script
########################################################################################


parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--mult_pep', type=str, help='multiple alignment in fasta format')
parser.add_argument('--mult_nuc', type=str, help='multiple alignment in fasta format')
parser.add_argument('--threshold', type=float, help='threshold used to remove unconfident alignment')
parser.add_argument('--guid_score', type=str, help='path to zorro')
parser.add_argument('--out', type=str, help='path to the output file')

args = parser.parse_args()
main(args.mult_pep, args.mult_nuc, args.guid_score, args.threshold, args.out)
