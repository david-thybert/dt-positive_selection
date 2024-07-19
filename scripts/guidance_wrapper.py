import argparse 
import subprocess
import os.path



"""
This script identify with zorro region of a multiple sequence alignment that is not confident and remove them.

command:

    python zorro_wrapper.py --zorro <path to zorro> --tree_cmd <path to Fastree> --mult <fmultiple alignment fasta> --threshold <0. - 10.> --out <filtered alignment> 

"""


def run_guidance(fasta_file:str, guidance_cmd:str, mafft_cmd:str, out_dir:str)->tuple:
    """
    Wrapper around zorro and return the score for the alignment in each position

    :param fasta_file: fast file of the sequence to align
    :param guidance_cmd: path  to guidance executable
    :param mafft_cmd: command to mafft (use absolute path for guidance)
    :return: the list of score corresponding at each column of the alignment
    """
    result = []
    command = f"{guidance_cmd} --seqFile {fasta_file} --msaProgram MAFFT --seqType aa --mafft {mafft_cmd}--outDir {out_dir}"
    out_process = subprocess.run(command, shell=True, capture_output=True)
    print(command)
    print(out_process.returncode)

    score_file = f"{out_dir}/MSA.MAFFT.Guidance2_col_col.scr"
    alignment_file = f"{out_dir}/MSA.MAFFT.aln.With_Names"

    if not os.path.exists(score_file):
        print(f"File {score_file} does not exist")
        raise FileExistsError(score_file)
    if not os.path.exists(alignment_file):
        print(f"File {alignment_file} does not exist")
        raise FileExistsError(alignment_file)

    return score_file, alignment_file


def main(fasta:str, guidance_cmd:str, mafft_cmd:str, out_dir:str)->None:
    """
    This is the main function of the script

    :param fasta: fasta file with the unaligned sequence fomr orthogroup
    :param guidance_cmd: path to guidance executable
    :param mafft_cmd: path to mafft executable
    :param out_dir: path where to store the out of guidance file
    :return: nothing
    """
    # run zorro to score alignment column confidence
    score_file, alignment_file = run_guidance(fasta, guidance_cmd, mafft_cmd, out_dir)
    gene_id = fasta.split("/")[-1].split(".pep")[0]
    os.rename(score_file, f"{gene_id}.ali_score")
    os.rename(alignment_file, f"{gene_id}.pep.ali")

                         
    

########################################################################################
########### Main script
########################################################################################


parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--fasta',type=str, help='fasta file of unaligned sequence fomr the orthogroup')
parser.add_argument('--guidance_cmd',type=str, help='guidance command')
parser.add_argument('--mafft_cmd', type=float, help='mafft command')
parser.add_argument('--out_dir', type=str, help='output directory')

args = parser.parse_args()
main(args.fasta, args.guidance_cmd, args.mafft_cmd, args.out_dir)

