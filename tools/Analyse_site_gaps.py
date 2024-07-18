
import argparse
import pandas as pd
from Bio import AlignIO



def load_position(possel:str, threshold)->dict:
    """

    """
    dico_pos = {}
    with open(possel) as file_handler:
        for line in file_handler:
            if "gene_id" in line:
                continue

            tab = line.split("\t")
            gene_id = tab[0]
            pos = tab[8]
            if pos == "":
                continue

            dico_pos[gene_id] = []
            tab_pos = pos.split("|")
            for pos_elem in tab_pos:
                pos = int(pos_elem.split(":")[0])
                score = float(pos_elem.split(":")[1])
                if score < threshold:
                    continue
                pos_id = f"{gene_id}-{pos}"
                dico_pos[gene_id].append([pos_id, pos])
    return dico_pos

def get_gap_profile(align_file)->list:
    """
    """
    print ("ici")
    with open(align_file) as file_handler:
        alignment = AlignIO.read(file_handler, "fasta")

    gap_list = []
    i = 0
    state = 0 # state no gap
    start_gap = 0
    end_gap = 0
    while i < len(alignment[0]):
        column = alignment[:, i]

        if "-" in column:
            if state == 0:# state no gap
                start_gap = i
                state = 1
            elif state == 1:# state gap
                pass
        else:
            if state == 1 :# state gap
                end_gap = i
                gap_list.append([start_gap, end_gap-1])
                state = 0
            elif state == 0:# state no gap
                pass
        i = i + 1
    return gap_list

def get_distance_gap_site(sites:list, gap_profile:list)->dict:
    """

    """
    result = {}
    for site in sites:
        pos_id = site[0]
        pos = site[1]
        min_dist = 10000000000
        for gap in gap_profile:
            if pos <= gap[1] and pos >=gap[0]:
                min_dist = 0
            else:
                ds = abs(pos - gap[0])
                de = abs(pos - gap[1])
                mind = de if de < ds else ds
                min_dist = mind if mind < min_dist else min_dist
        result[pos_id] = min_dist
    return result
def main(possel:str, dir_align:str,thr:float, out:str)->None:
    """
    
    """
    dico_pos = load_position(possel, thr)
    print(dico_pos)
    result = {}
    for ortho_gene, positions in dico_pos.items():
        align_file = f"{dir_align}/{ortho_gene}.pep.fasta.ali"
       #align_file  = "/Users/thybed/Documents/Workspace/Analysis/ENSADMT00000000906.1.pep.fasta.ali.fa"
        gap_profile = get_gap_profile(align_file)
        print(gap_profile)
        dico_dist = get_distance_gap_site(positions, gap_profile)
        result.update(dico_dist)
        print(result)

    with open(out, "w") as file_handler:
        for pos_id, dist in result.items():
            file_handler.write(f"{pos_id}\t{dist}\n")



parser = argparse.ArgumentParser(description='Script analysing the site under positive selction and their distance to gaps')
parser.add_argument('--possel',type=str, help='Possel output file with significant hits')
parser.add_argument('--dir_align', type=str, help='directory wher pep alignment are located')
parser.add_argument('--thr', type=float, default=0.95, help='threshold for analysinf position')
parser.add_argument('--out', type=str, help='outfile')


args = parser.parse_args()
main(args.possel, args.dir_align, args.thr, args.out)