import argparse 
from Bio import SeqIO


def rename_id_pep(pep_file:str)->list:
    results = []
    dico_genes = {}
    seq_records = SeqIO.parse(pep_file, "fasta")
    for record in seq_records:
        lst_item = record.description.split()
        for item in lst_item:
            if "transcript:" in item:
                id_tr = item.split(":")[-1]
            if "gene:" in item:
                id_gene = item.split(":")[-1]
        if not id_gene in dico_genes:
            dico_genes[id_gene] = []
        record.description = id_tr
        record.id = id_tr
        results.append(record)
    return results


def rename_id_tr(cds_file:str, peps:list)->list:
    results = []
    dico_tr = {}
    seq_records = SeqIO.parse(cds_file, "fasta")
    for record in seq_records:
        record.decription = record.id
        results.append(record)
    return results


def main(cds_file:str, pep_file:str, out_cds:str, out_pep:str)-> None:
    peps = rename_id_pep(pep_file)
    print(len(peps))
    trs = rename_id_tr(cds_file, peps)
    print(len(trs))
    with open(out_cds, "w") as output_handle:
        SeqIO.write(trs, output_handle, "fasta")
    with open(out_pep, "w") as output_handle:
        SeqIO.write(peps, output_handle, "fasta")

###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--cds',type=str, help='file describing the orhtology relationship between genes ')
parser.add_argument('--pep', type=str, help='path to the directory contianing the CDS files')
parser.add_argument('--out_cds', type=str, help='path to the directory contianing the peptide files')
parser.add_argument('--out_pep', type=str, help='path to the outdirectory')

args = parser.parse_args()
main(args.cds, args.pep, args.out_cds, args.out_pep)

