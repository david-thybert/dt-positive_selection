import sys
import argparse

def main(ortho_file:str, min_species:int, required_species:str, out_file:str )-> None:
    
    dico_ortho = {}
    species = []

    req_species = required_species.split()

    with open(ortho_file) as file_handler:
        for line in file_handler:
            if "Orthogroup" in line:
                species = line.split("\t")
                continue
            ortho = line.split("\t")
            i = 1
            dico_ortho[ortho[0]] = [""]*len(species)
            while i < len(ortho):
                tab_elem = ortho[i].split(",")
                if len(tab_elem) == 1:
                    dico_ortho[ortho[0]][i-1] = tab_elem[0]
                i = i + 1

    with open(out_file, 'w') as outfile_handler:
        outfile_handler.write("\t".join(species[1:]) + "\n")
        for ortho, lst_ortho  in dico_ortho.items():
            nb = 0
            for gene in lst_ortho: 
                if gene != "":
                    nb = nb + 1
            if nb >= min_species:
                outfile_handler.write("\t".join(lst_ortho))





###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--ortho',type=str, help='file describing the orhtology relationship between genes ')
parser.add_argument('--min_species', type=int, help='path to the directory contianing the CDS files')
parser.add_argument('--req_sp', type=str, help='path to the directory contianing the peptide files')
parser.add_argument('--out', type=str, help='path to the outdirectory')

args = parser.parse_args()
main(args.ortho, args.min_species, args.req_sp, args.out)
