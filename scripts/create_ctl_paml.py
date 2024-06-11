import argparse

"""
This script create a ctl file describing the model to run using codeml. The model implemented are the null model for branch site test (fixing w to 1)
and the alternative branch site test where w is estimated from the data.

command:

    python create_ctl_paml.py --type <null|alt> --alignment <alignemnt file .phy> --tree <tree file (with tag)> --out_PAML <path to cdt file> --out_file <path to ctl file>

"""

def create_empty_ctl_dic()->dict:
    """
    This function create an empty dict will all filed required by ctl file

    :return : dictionary with inital values
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
    Main function of th script
    
    :param  alignment: path to lignment file in phylip format
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

parser = argparse.ArgumentParser(description='This script create aconfiguration file for codeml model')
parser.add_argument('--alignment',type=str, help='path to multiple sequence alignment in phylip format')
parser.add_argument('--tree',type=str, help='path tot he tree in nwk format *(tree need to be tagged)') 
parser.add_argument('--out_PAML',type=str, help='path to the future cdt file genereated by codeml')
parser.add_argument('--out_file',type=str, help='path to the output ctl file ')
parser.add_argument('--type',type=str, help='null|alt whther the script generate a confidguraiton for the null or alternative model ')

args = parser.parse_args()
main(args.alignment, args.tree, args.type, args.out_PAML, args.out_file)