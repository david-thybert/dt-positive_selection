import sys
import os
import argparse 


"""
This script split an orthogroup in batches of a size gevin by the user
The script generate a set of files corresponding to a batch with the same format as the outhogroup file

command:

    python batch_orthogroups.py --ortho < orthgroup file> --size <integer> --out_prefix <prefix for ouput file>

"""

def load_ortho_file(orthogroups:str)->tuple:
    """
    This function load the orthgorup file and separate the header fomr the content

    :param orthogroups: the path to the orthogroup file
    :return: a tuple composed of the list of lines split into a list and the head split into a list
    """
    
    species = []
    lst_orthologues = []

    with open(orthogroups) as file_handler:
        lst_lines = file_handler.readlines()
        species = lst_lines[0].replace("\n","").split()
        for line in lst_lines[1:]:
            if line.strip() == "":
                continue
            orthologues = line.replace("\n","").split("\t")
            lst_orthologues.append(orthologues)
    
    return (lst_orthologues, species)

def create_batches(lst_orthologues:list, species:list, size:int)->list:
    """
    split the list of orthologues into batches

    :param lst_orthologues: list of othologues arrays
    :param species: species header
    :param size: size of the batch
    :return : list of batches compaosed of a species header and number of orthologues of size size 
    """

    batches = []
    batch = []
    i = 0
    while i < len(lst_orthologues):
        
        if len(batch) == 0:
            batch.append(species) # add header
        
        batch.append(lst_orthologues[i])
        
        if len(batch) == size + 1 : # size + 1 because counting the header in the list
            batches.append(batch)
            batch = []
        
        i = i + 1

    if batch != []:
        batches.append(batch)

    return batches

def main(orthogroups:str, size:int, out_prefix:str)->None:
    """
    Main function of the script

    :param orthogroups: path to the orthgourp file
    :param size: size of the batch
    :param out: outdirectory
    """
    (lst_orthologues, species) = load_ortho_file(orthogroups)

    batches = create_batches(lst_orthologues, species, size)

    index = 0
    for batch in batches:
        path_file = f"{out_prefix}_{index}.obatch"
        with open(path_file, 'w') as file_handler:
            for elem in batch:
                file_handler.write("\t".join(elem) + "\n")
        index = index + 1
    
########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script spliting the orthogroup files in batches')
parser.add_argument('--ortho',type=str, help='orthogroup file')
parser.add_argument('--size', type=int, help='batch size')
parser.add_argument('--out_prefix', type=str, help='outfile prefix')

args = parser.parse_args()
main(args.ortho, args.size, args.out_prefix)

