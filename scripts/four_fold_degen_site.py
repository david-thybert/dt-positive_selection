import typing
import argparse 

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
from Bio import AlignIO

def four_fold_degen(triplet: str) -> str:
    '''
    This method return the last nucleotide of a fou fold degenrate site codon.
    if the codon is not four fold degenratethen return empty string.
    "CGU, CGC, CGA, CGG" --> Arg
    "GCU, GCC, GCA, GCG" --> Ala
    "GGU, GGC, GGA, GGG" --> Gly
    "CUU, CUC, CUA, CUG" --> Leu
    "UCU, UCC, UCA, UCG" --> Ser
    "ACU, ACC, ACA, ACG" --> Thr
    "GUU, GUC, GUA, GUG" --> Val 

    :param triplet: codon to be evaluated
    :return: returns a nucleotide if the codon is four fold degenrated , an empty string otherwise
    '''
    if triplet[-1] == '-':
        return ""
    if triplet[0:2] == "CG":# Arg
        return triplet[-1]
    elif triplet[0:2] == "GC":# Ala
        return triplet[-1]
    elif triplet[0:2] == "GG":# Gly
        return triplet[-1]
    elif  triplet[0:2] == "CT":# Leu
        return triplet[-1]
    elif triplet[0:2] == "TC":# Ser
        return triplet[-1]
    elif triplet[0:2] == "AC":# thr
        return triplet[-1]
    elif triplet[0:2]=="GT":# Val
        return triplet[-1]
    return ""


def _to_alignment_object(matrix: list, ori_align:Alignment) -> Alignment:
    '''
    This method build a multiple alignment object from the four fold degenated site matrix 
    and keep the same sequecne id that in the original multiple alignment

    :param matrix: matrix of four fold degenrated site
    :ori_align: original multiple sequence alignment
    :return: multiple sequence alignment object of four fold degenerated sites.
    '''
    sequences = []
    i = 0
    while i  < len(matrix[0]):
        strSeq = ""
        for column in matrix:
            strSeq = strSeq + column[i]
        sequence = SeqRecord(Seq(strSeq), id=ori_align[i].id)
        sequences.append(sequence)
        i = i + 1
    coordinates = Alignment.infer_coordinates(sequences)
    alignment = Alignment(sequences, coordinates)
    return alignment


def get_4_fold_dege_align(alignment: Alignment) -> Alignment:
    '''
    Retrive all four fold degenrate site that ahs the two first 
    nucleotid of the codon conserved across the species analysed

    :param alignment: multiple sequence alignment of the cDNA at the nucleotide level
    :return: multiple alignment objecit of alny the four fold degenerated sites
    '''
    index_codon = 0
    deg_align = []
    while index_codon < len(alignment[0]):
        degen = True
        conserved = True
        prevCodon = ""
        deg_column = []
        for align_seq in alignment:
            if prevCodon != "" and prevCodon[0:-1] != align_seq[index_codon:index_codon+2]: # we consider only conserved 4 fold degenerate codon on the two first nuc.
                conserved = False
                break
            nuc = four_fold_degen(str(align_seq[index_codon:index_codon+3].seq))
            deg_column.append(nuc)
            if nuc == "": # if not degenerated
                degen = False
                break
        if degen and conserved:
            deg_align.append(deg_column)
        index_codon = index_codon + 3   
    
    four_fold_align = _to_alignment_object(deg_align, alignment)
    print(four_fold_align)
    return four_fold_align


def main(nuc_align: str, outfile: str) -> None:
    '''
    Main fucntiona of the script

    :param nuc_align: path tot he multiple sequence alignment
    :param outfile: path tothe out file storing the four fold degenrated multiple alignment.
    '''
    alignment = AlignIO.read(open(nuc_align), "fasta")
    four_fold_align = get_4_fold_dege_align(alignment)
    with open(outfile, "w") as file_hanlder:
        file_hanlder.write(format(four_fold_align, "fasta"))
    


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='The script check that the nucleotide unaligned file correspsond to the aligned peptide file')
parser.add_argument('--nucal',type=str, help='peptide alignemnt file')
parser.add_argument('--out', type=str, help='path to the outdirectory')
args = parser.parse_args()

main(args.nucal, args.out)