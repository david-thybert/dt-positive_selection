import argparse 
from Bio import SeqIO



def _get_longest(lst_records:list)-> object:
    record  = lst_records[0]
    for rec in lst_records:
        if len(record) < len(rec):
            record = rec
    return record


def select_pep(pep_file:str)->list:
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
        dico_genes[id_gene].append(record)

    for gene, peps in dico_genes.items():
        longest_pep = _get_longest(peps)
        results.append(longest_pep)
    return results


def select_transcript(cds_file:str, peps:list)->list:
    results = []
    dico_tr = {}
    seq_records = SeqIO.parse(cds_file, "fasta")
    for record in seq_records:
        record.decription = record.id
        dico_tr[record.id] = record

    for pep in peps:
        results.append(dico_tr[pep.id])
    return results


def main(cds_file:str, pep_file:str, out_cds:str, out_pep:str)-> None:
    peps = select_pep(pep_file)
    print(len(peps))
    trs = select_transcript(cds_file, peps)
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

