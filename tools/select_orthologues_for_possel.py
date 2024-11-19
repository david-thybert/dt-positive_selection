import sys
import argparse
import json
from pprint import pprint
import os

def load_conf(conf:str)->dict:
    """
    """
    result = {}
    with open(conf) as f:
        result = json.load(f)
    return result

def get_species_file_mapping(dir:str)->dict:
    """
    """
    result = {}
    files = os.listdir(dir)
    for file in files:
        species = file.split("__v__")[-1].split(".")[0]
        result[species] = f"{dir}/{file}"
    return result

def _get_empty_gene_dict(spe_to_consider):
    """
    """
    result = {}
    for species in spe_to_consider: result[species]=[]
    return result

def create_gene_spe_matrix(spe_to_consider:list, spe_file_map:dict)->dict:
    """
    """
    result = {}
    for species in spe_to_consider:
        with open(spe_file_map[species]) as file_hand:
            for line in file_hand:
                if "Orthogroup" in line:
                    continue
                tab  = line.strip().split("\t")
                genes_ref = tab[1].split(",") 
                if len(genes_ref) > 1: # we need only 1-to-1 orthologues 
                    continue
                if not genes_ref[0] in result:
                    result[genes_ref[0]] = _get_empty_gene_dict(spe_to_consider)
                result[genes_ref[0]] = tab[2].replace(" ","").split(",") 
    return result

def filter_genes_from_constraint(gene_matrix:dict, conf:dict)->dict:
    """
    """
    result = {}
    # for each gene removing duplicate in each species
    for gene, spe_dict in gene_matrix.items():
        for species, othologues in spe_dict.items():
            if len(othologues)> 1:
                gene_matrix[gene][species]=[]
    
    # select only gene that complies to the constraint
    for gene, spe_dict in gene_matrix.items():
        
        species_with_gene = []
        for species, othologues in spe_dict.items():
            if len(othologues)==1:
                species_with_gene.append(species)

        # constraint of size
        if len(species_with_gene) < conf["min_species"]:
            continue

        # constraint on compulsory
        tag_present = True
        for species in conf["compulsory"]:
            if not species in species_with_gene:
                tag_present = False
                break
        if not tag_present: continue

        result[gene] = gene_matrix[gene]
    return result


def save_matrix(matrix:dict, conf:dict, outfile:str)-> None:
    
    # we want the reference species to out of this list
    lst_species = [species for species in conf["species"] if species != conf["ref"] ]
    # we want the ref species to be the first species in the header
    lst_species_header = [conf["ref"]] + lst_species
    header = "\t".join(lst_species_header)

    with open(outfile, "w") as out_handler:
        out_handler.write(f"{header}\n")
        for gene, spe_dict in matrix.items():
            gene_tab = [gene] + [spe_dict[species][0] if len(spe_dict[species])==1  else "" for species in conf["species"] ]
            out_handler.write(f"{"\t".join(gene_tab)}\n")


def main(dir_ref:str, conf:str, out:str)->None:
    """
    """
    # load the configuraiton file
    conf = load_conf(conf)

    # make a species-file mapping
    spe_file_map = get_species_file_mapping(dir_ref)
    
    # create a gene species matric based in configuraiton oinfo and the specoes-file mapping
    gene_spe_matrix = create_gene_spe_matrix(conf["species"], spe_file_map)

    # remove the gene form the matrix that does not satisfy the constraint
    filtered_spe_matrix = filter_genes_from_constraint(gene_spe_matrix, conf)
    
    # save the matrix
    save_matrix(filtered_spe_matrix, out)



###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--dir',type=str, help='file describing the orhtology relationship between genes ')
parser.add_argument('--conf', type=str, help='path to the directory contianing the CDS files')
parser.add_argument('--out', type=str, help='path to the outdirectory')

args = parser.parse_args()
main(args.dir, args.conf,args.out)

