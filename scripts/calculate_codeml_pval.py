
import argparse 
from scipy.stats.distributions import chi2

def parse_codeml_output(ctd_file:str)-> dict:
    """
    """
    result = {}
    with open(ctd_file) as file_hanlder:
        lines = file_hanlder.readlines()

    for line in lines:
        if "lnL(" in line: # likelihood line
            np = int(line.split(")")[0].split(":")[-1])
            result["np"] = np
            log_lh = float(line.split(":")[-1].strip().split()[0])
            result["log_lh"] = log_lh
        if  "background w" in line:
            back_w_0 = float(line.split()[2])
            result["back_w_0"] = back_w_0
            back_w_1 = float(line.split()[3])
            result["back_w_1"] = back_w_1
            back_w_2a = float(line.split()[4])
            result["back_w_2a"] = back_w_2a
            back_w_2b = float(line.split()[5])
            result["back_w_2b"] = back_w_2b
        if  "foreground w " in line:
            fore_w_0 = float(line.split()[2])
            result["fore_w_0"] = fore_w_0
            fore_w_1 = float(line.split()[3])
            result["fore_w_1"] = fore_w_1
            fore_w_2a = float(line.split()[4])
            result["fore_w_2a"] = fore_w_2a
            fore_w_2b = float(line.split()[5])
            result["fore_w_2b"] = fore_w_2b
    return result

def main(alt_codeml_file:str, null_codeml_file:str, out_file:str)->None:
    
    # parse results from file
    alt_codeml = parse_codeml_output(alt_codeml_file)
    null_codeml = parse_codeml_output(null_codeml_file)

    # calculate pval that alternative model is more likely
    df  = alt_codeml["np"] - null_codeml["np"]
    delta_lrt = 2*(alt_codeml["log_lh"] - null_codeml["log_lh"]) 
    pval =  chi2.sf(delta_lrt, df)

    # get the gene id fomr the file name
    gene_id = alt_codeml_file.split("/")[-1].split(".alt")[0]
    
    # save output
    header = "gene_id\tlrt\tpval\talt_w_0\talt_w_1\talt_w_2a\talt_w_2b\tnull_w_0\tnull_w_1\tnull_w_2a\tnull_w_2b"
    line = f"{gene_id}\t{delta_lrt}\t{pval}\t{alt_codeml['fore_w_0']}\t{alt_codeml['fore_w_1']}\t{alt_codeml['fore_w_2a']}\t{alt_codeml['fore_w_2b']}"
    line = line + f"\t{null_codeml['fore_w_0']}\t{null_codeml['fore_w_1']}\t{null_codeml['fore_w_2a']}\t{null_codeml['fore_w_2b']}"

    with open(out_file, "w") as file_handler:
        file_handler.write(f"{header}\n")
        file_handler.write(f"{line}\n")

########################################################################################
########### Main script
########################################################################################


parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--alt_codeml',type=str, help='multiple alignment in fasta format')
parser.add_argument('--null_codeml', type=str, help='threshold used to remove unconfident alignment')
parser.add_argument('--out', type=str, help='path to the output file')

args = parser.parse_args()

main(args.alt_codeml, args.null_codeml,args.out)