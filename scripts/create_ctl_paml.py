import argparse

def create_empty_ctl_dic()->dict:
    """
    """
    ctl = {"seqfile":"", "treefile":"", "outfile":"",
           "noisy":9, "verbose":1,"runmode":0,
           "seqtype":1,"CodonFreq":2,"clock":0, 
           "aaDist":0, "model":2, "NSsites":2,
           "icode":0, "Mgene":0, "fix_kappa":0,
           "kappa":2, "fix_omega":1, "omega":1,
           "getSE":0, "RateAncestor":0, "Small_Diff":0.45e-6,
           "cleandata":1, "fix_blength":0
           }
    return ctl

def main(alignment:str, tree:str, type:str, out_PAML:str, out_file:str)-> None:
    """
    """
    ctl = create_empty_ctl_dic()

    ctl["seqfile"] = alignment
    ctl["treefile"] = tree
    ctl["outfile"] = out_PAML

    if type == "null":
        ctl["fix_omega"] = 1
    elif type == "alt":
        ctl["fix_omega"] = 0
    else:
        print(f"--type parameter is wrong : {type} not recognised only 'null' or 'alt' allowed")
        sys.exit(1)

    with open(out_file, "w") as outfile_handler:
        for param,val in ctl.items():
            outfile_handler.write(f"{param} = {val} \n")
    


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='The script check that the nucleotide unaligned file correspsond to the aligned peptide file')
parser.add_argument('--alignment',type=str, help='list of degenerated files')
parser.add_argument('--tree',type=str, help='base directory for the files') 
parser.add_argument('--out_PAML',type=str, help='base directory for the files')
parser.add_argument('--out_file',type=str, help='base directory for the files')
parser.add_argument('--type',type=str, help='null|alt')
args = parser.parse_args()

main(args.alignment, args.tree, args.type, args.out_PAML, args.out_file)