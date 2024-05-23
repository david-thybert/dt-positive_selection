import typing
import sys
import os
import argparse 
from Bio import SeqIO

def check_correspondance(four_folds: list, trees: list)-> list:
    """
    Look if each four fold degenrated alignemnt have a 
    corresponding tree file

    :param four_folds: list fo four fold degenerated alignemnt files
    :param trees:list of tree files
    :return:return a list of tuple with corespsonding four fold degen alignemnt and tree.[(four_fold, tree),...]
    """
    results=[]
    for ff_elem in four_folds:
        for tree_elem in trees:
            if ff_elem.split(".")[0] == tree_elem.split(".")[0]:
                results.append((ff_elem,tree_elem))
                break
    return results   

def check_tree_ok(trees:list)-> (list):
    """
    to implment
    """
    pass

def check_ff_ok(four_folds)->(list):
    """
    to implement
    """
    pass


def main(ff_align_lst_file: str, tree_lst_file: str, out_base_dir: str, out_file:str)-> None:
    """
    Main funciton fo the script, check that the two sets of list are without error and corresponds

    :param ff_align_lst_file: list of file name for the four fold degenrate site alignemnt 
    :param tree_lst_file: list fo file name for the tree files
    :out_base_dir: base directory of the outfile
    :out_file: name of the outfile
    :return: returns nothing
    """
    
    ff_files = ff_align_lst_file.split(" ")
    tree_files = tree_lst_file.split(" ")

    pair_files = check_correspondance(ff_files, tree_files)
    with open(out_file, "w") as file_handler:
        print("wite files\n")
        for ff_deg_tree_pair in pair_files:
            file_handler.write(f"{ff_deg_tree_pair[0]}\t{ff_deg_tree_pair[1]}\n")



parser = argparse.ArgumentParser(description='The script check that the nucleotide unaligned file correspsond to the aligned peptide file')
parser.add_argument('--ff_deg',type=str, help='list of degenerated files')
parser.add_argument('--base_dir',type=str, help='base directory for the files') 
parser.add_argument('--trees',type=str, help='list of trees')
parser.add_argument('--out', type=str, help='outfiles')
args = parser.parse_args()

main(args.ff_deg, args.trees, args.base_dir, args.out)
