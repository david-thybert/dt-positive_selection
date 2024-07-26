import argparse
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
from Bio import AlignIO
import pandas as pd

"""
This script filters alignemnt based on list of column quality alignment score and gaps form the alignemnt.

command:

    python filter_alignment.py --mult_pep <pep multiple alignment> --mult_nuc <nuc multiple alignment> 
                               --guid_score <guidance score> --threshold <0 - 1> 
                               --neighbors <0,1,...,n> --out <filtered alignment> 

"""

def load_align(mult_fasta: str) -> object:
    """
    Load multiple peptide sequence alignment

    :param mult_fasta:  multiple alignment file in fasta format
    :return: Alignment object
    """
    alignment = None
    with open(mult_fasta) as file_handler:
        alignment = AlignIO.read(file_handler, "fasta")
    return alignment

def get_gapped_regions(align: object, neighbor:int)->list:
    """
    Return a description of gapped and non-gapped part of the alignment,
    It will consider as gapped the +/- neighbor around a gap

    :param align: the alignment object
    :param neighbor: the number of neighbor to identify as gapped in the alignment.
    :return: list of the size of the alignment masking with 0 the position with a gapped
    """
    i = 0
    result = [1] * len(align[0])
    state = 0
    start = -1
    end = 0
    while i < len(result):
        print(start, end)
        print(align[:,i])
        if "-" in align[:,i]:
        #    print("gap")
            if state == 0:
                start = i - neighbor if (i-neighbor) >=0 else 0
                state = 1
            if state == 1:
                end = i
        else:
            if state == 1:
                end = i 
                j = start
                while j < end:
                    result[j] = 0
                    j=j+1
                state = 0
                start = -1
                end =0
        i = i +1
    if state == 1:
        j = start
        while j <= end:
            result[j] = 0
            j = j + 1
    return result



def get_confident_regions(scored_pos: list, gapped_pos: list, threshold: float) -> list:
    """
    Get a list fo region of confident alignment from column alignment score and gapps
    This function works at the peptide level

    :param scored_pos: list fo confidence score per position in the alignment
    :param gapped_pos: list of gapped/non-gapped masking per position in the alignment
    :param threshold: threshold define a confident/non-confident alignment
    :return: list of regions [start, end] that are confidently well aligned and without gaps
    """
    i = 0
    state = 0  # state 0 : belongs to a confident aligned position
    # state 1 : belong to a non-confident aligned position
    start = 0
    end = -1
    confident = []
    while i < len(scored_pos):
        if scored_pos[i] < threshold or gapped_pos[i] == 0:
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
    new alignment object

    :param alignment: alignment object
    :param confidence_pos: list of confident regions in the format of [[start:end],...]
    :return: concatenation of the confident regions in the order in an alignment object
    """
    filtered_alignment = None

    if confidence_pos == []:
        return filtered_alignment

    for interval in confidence_pos:
        start = interval[0]
        end = interval[-1]
        if filtered_alignment is None:
            filtered_alignment = alignment[:, start:end]
        else:
            filtered_alignment = filtered_alignment + alignment[:, start:end]

    return filtered_alignment


def filter_align_nuc(alignment: object, confidence_pos: list) -> object:
    """
    The function keep the high confident column of the alignment in a
    new nucleotide alignment object

    :param alignment: alignment object
    :param confidence_pos: list of confident regions in the format of [[start:end],...] in aminoacid coordinate
    :return: concatenation of the confident regions in the order in an alignment object
    """
    filtered_alignment = None

    if confidence_pos == []:
        return filtered_alignment

    for interval in confidence_pos:
        start = interval[0] * 3
        end = interval[-1] * 3
        if filtered_alignment is None:
            filtered_alignment = alignment[:, start:end]
        else:
            filtered_alignment = filtered_alignment + alignment[:, start:end]

    return filtered_alignment

def load_score(score_file: str)-> list:
    """
    load score from alignment column quality file

    :param score_file: file describing the alignment
    :return: list with a score for each alignment position
    """
    result = []
    with open(score_file) as file_hand:
        for line in file_hand:
            if "#" in line:
                continue
            tab = line.split()
            result.append(float(tab[1]))
    return result


def main(mult_pep_fasta: str, mult_nuc_fasta: str, score_file: str, threshold: float, neighbors: int, out_file: str) -> None:
    """
    This is the main function of the script

    :param mult_pep_fasta: peptide multiple alignment in fasta format
    :param mult_nuc_fasta: nucleotide multiple alignment in fasta format
    :param score_file: command for tree tool
    :param threshold: threshold for filtering no confident alignment
    :param neighbors: the +/- n neighbor around a gap that need to be filtered.
    :param out_file: out_File storing the filtered alignment
    :return: nothing
    """
    # identify highly confident position based on threshold
    scored_pos = load_score(score_file)
    alignment_pep = load_align(mult_pep_fasta)
    gapped_pos = get_gapped_regions(alignment_pep, neighbors)
    confidence_pos = get_confident_regions(scored_pos, gapped_pos, threshold)
    print(alignment_pep)
    print(gapped_pos)
    if confidence_pos == []:
        return

    # keep only highly confident regions from alignment
    alignment_pep_filtered = filter_align_pep(alignment_pep, confidence_pos)
    alignment_nuc = load_align(mult_nuc_fasta)
    alignment_nuc_filtered = filter_align_nuc(alignment_nuc, confidence_pos)
    print(alignment_pep_filtered)
    print(alignment_nuc_filtered)
    # save the filtered alignments
    with open(f"{out_file}.pep.filt", "w") as file_hanlder:
        file_hanlder.write(format(alignment_pep_filtered, "fasta"))
    with open(f"{out_file}.nuc.filt", "w") as file_hanlder:
        file_hanlder.write(format(alignment_nuc_filtered, "fasta"))

    # save the regions files
    regions_file = f"{out_file}.confident_reg"
    with open(regions_file, "w") as file_handler:
        for pos in confidence_pos:
            file_handler.write(f"{pos[0]}\t{pos[-1]}\n")


########################################################################################
###########  Main script
########################################################################################


parser = argparse.ArgumentParser(description='Script filtering low quality alignment column and gapped regions')
parser.add_argument('--mult_pep', type=str, help='multiple alignment in fasta format')
parser.add_argument('--mult_nuc', type=str, help='multiple alignment in fasta format')
parser.add_argument('--threshold', type=float, help='threshold used to remove unconfident alignment')
parser.add_argument('--guid_score', type=str, help='path to zorro')
parser.add_argument('--neighbors', type=int, default=0, help='number of n neighbor to be removed with gaps')
parser.add_argument('--out', type=str, help='path to the output file')

args = parser.parse_args()
main(args.mult_pep, args.mult_nuc, args.guid_score, args.threshold, args.neighbors, args.out)
