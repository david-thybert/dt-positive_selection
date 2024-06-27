
import typing
import sys
import os
import argparse 
from Bio import SeqIO


"""
This script format the fasta file from a fasta file per species will all the proteiome or transcriptiome 
to a fasta file per orthogroup for proteome and transcriptome. To do this it uses an orhtology matrix 
describing 121 ortholgous shared by the species in the patrix. 

command:

    python format_fasta.py --ortho <ortholog matrix> --nuc <nuc direcotry> --pep <pep directory> --out <out directory>

"""


def _getTranscripr_id(description: str) -> str:
    """
    Read an Ensembl peptide fasta header and return the transcript ID

    :param description: the header of the fasta (description field of a biopython record)
    :return: returns the Ensembl transcript id 
    """
    lst_descripton  = description.split()
    for description in lst_descripton:
        if "transcript:" in description:
            return description.split(":")[-1]
    return ""


def readOrthologyFile(orthologues: str) -> dict:
    """
    Read an orthology matrix file and return dictionary obejct representing the matrix.
    In the matrix file, each column is a species and each row is a orthogroup. Each cell
    contains the id of the cannonical trasncript for the orthologous gene for a given 
    species and a given orthologous group. The transcript id of the first species will 
    be used to name the orthogroup

    :param orthologues: the path tot he ortholog matrix file
    :return: return a dictionary representing the homologius matrix {horthologous group: {species:[transcriptid,'','']}}.
    
    """
    result = {}
    with open(orthologues) as file_handler:
        species_list = next(file_handler).strip().split("\t")
        for line in file_handler:
            lst_gene = line.strip().split("\t")
            if lst_gene == ['']:
                break
            result[lst_gene[0]] = {}
            i = 0
            while i < len(species_list):
                result[lst_gene[0]][species_list[i]] = [lst_gene[i],"",""]
                i = i + 1
    return result


def createPepFile(orthogroup_path: str, orthogroup: str, peptide: str, orthodata: dict) -> list :
    '''
    This method group all peptide sequences from each orthogroup into one file per orthogroup. 
    {horthologous group: {species:[transcriptid,'','']}}.

    :param orthogroup_path: basal path where the orthogroup files will be stored
    :param orthogroup: name of the orthogroup
    :param peptide: basal path where all the peptide file are saved for each species
    :param orthodata: dictionary containing the orthology matrix
    :return: returns nothing
    '''
    fasta_list = []
    for species, genes in orthodata.items():
        species_sequences = f"{peptide}/{species}.fa"
        seq_records = SeqIO.parse(species_sequences, "fasta")
        for record in seq_records:
            #if _getTranscripr_id(record.description) == genes[0]:
            if record.id == genes[0]:
                description  = f"species:{species} transcript:{genes[0]} ortho:{orthogroup}"
                record.id = f"{genes[0]}|{species}"
                record.description = ""
                fasta_list.append(record)
    return fasta_list


def createNucFile(orthogroup_path: str, orthogroup: str, nucleotide: str, orthodata: dict) -> list :
    """
    This method group all nucleotide sequences from each orthogroup into one file per orthogroup. 

    :param orthogroup_path: basal path where the orthogroup files will be stored
    :param orthogroup: name of the orthogroup
    :param nucleotide: basal path where all the nucleotide files are saved for each species
    :param orthodata: dictionary containing the orthology matrix
    :return: returns nothing
    """
    fasta_list = []
    for species, genes in orthodata.items():
        species_sequences = f"{nucleotide}/{species}.fa"
        print(species_sequences)
        record_dict = SeqIO.to_dict(SeqIO.parse(species_sequences, "fasta"))
        try:
            record = record_dict[genes[0]]
            description  = f"species:{species} transcript:{genes[0]} ortho:{orthogroup}"
            record.id = f"{genes[0]}|{species}"
            record.description = ""
            fasta_list.append(record)
        except KeyError as err:
            print(f"{genes[0]} absent from the fasta file : {err}")
    return fasta_list

def check_consistancy(lst_pep_fasta, lst_nuc_fasta) -> bool:
    """
     This method check that CDNA and pep sequence corresponde

     :param lst_pep_fasta: list fo peptide sequence
     :param lst_nuc_fasta: list of nucleotide sequences
     :return: True if all nuc sequence correpsont to the peptide equivalent False otherwise 
    """
    return True
    for nuc in lst_nuc_fasta:
        found = False
        if len(nuc) %3 != 0:
            return False
        
        for pep in lst_pep_fasta:
            if nuc.id == pep.id:
                found = True
                nuc_translate = nuc.seq.translate(to_stop=True)
                #print(str(nuc_translate))
                #print(str(pep.seq))
                if str(nuc_translate) == str(pep.seq):
                    break
                return False
        if found == False: 
            return False
    return True 

def main(orthologues: str, nucleotide: str, peptides: str, outBase: str) -> None:
    """
    Main function of the script. It reads the orhtology matrix and then use 
    the info from the orthology matrix to group pep and nuc sequenes in a
    coprrepsonding pep and nuc orhtogroup file

    :param orthologues: path to the orthologous matrix
    :param nucleotide: path tot he nucleotid file
    :param peptides: path to the peptide file
    :param outBase: path to the out directory
    :return: returns nothing
    """
   
    try :
        os.mkdir(outBase)       
    except FileExistsError:
        print(f"{outBase} Directory already exists.")
    except OSError as err:
        print(f"Error creating directory: {err}")
        sys.exit(1)
    
    ortho_data  = readOrthologyFile(orthologues)
    for ortho, lst_genes in ortho_data.items():
        fasta_pep = createPepFile(outBase, ortho, peptides, lst_genes)
        fasta_nuc = createNucFile(outBase, ortho, nucleotide, lst_genes)
        if check_consistancy(fasta_pep, fasta_nuc):
            nuc_file = f"{outBase}/{ortho}.nuc.fasta"
            pep_file = f"{outBase}/{ortho}.pep.fasta"
            with open(nuc_file, "w") as output_handle:
                SeqIO.write(fasta_nuc, output_handle, "fasta")
            with open(pep_file, "w") as output_handle:
                SeqIO.write(fasta_pep, output_handle, "fasta")
        else:
            print(f"issue in the orthogroup {ortho}: the nuceotide file does not correpsondant to pep file")    


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--ortho',type=str, help='file describing the orhtology relationship between genes ')
parser.add_argument('--nuc', type=str, help='path to the directory contianing the CDS files')
parser.add_argument('--pep', type=str, help='path to the directory contianing the peptide files')
parser.add_argument('--out', type=str, help='path to the outdirectory')

args = parser.parse_args()
main(args.ortho, args.nuc, args.pep, args.out)


