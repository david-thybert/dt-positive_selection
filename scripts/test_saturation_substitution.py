
import math
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind
from Bio import AlignIO
import argparse 

"""
This script test if the overall sites within an alignemnt show evidence that the sequecne is not saturated
It compare the distribution of obserevd entropy accross sites vs the expected entropy under saturation.
The expected entropy under saturation is cacluated based on the overal frequency of each nucleotide class accorss 
all infromative sites (have at least two different nuc class represented in the site and no gaps)

command:

    python test_saturation_substitution.py -mult <alignment.fasta> --out <outfile.tsv>

"""


def load_align(mult_fasta:str)->object:
    """ Load multiple peptide sequence alignment

    :param mult_fasta:  multiple alignement file in fasta format
    :return: Alignment object
    """
    alignment = None
    with open(mult_fasta) as file_handler:
        alignment = AlignIO.read(file_handler, "fasta")
    return alignment

def calculate_entropy(nb_seq:int, p:dict)-> float:
    """ This function calculate the entropy based on 
    https://doi.org/10.1093/sysbio/syab075.

    :param nb_seq: number of sequence in the alignment
    :param p: dictionary given the frequency of apparition 
              for each classe {A,T,G,C}.
    :return: the entropy value for a given class frequecny distribution 
    """
    nb_classes = len(p)
    
    # Step 1: Calculate log(n!)
    log_n_fact = math.log(math.factorial(nb_seq))
    
    # Step 2: Calculate -n * sum(p_i * log(p_i))
    entropy_term = -nb_seq * sum(p_i * math.log(p_i) if p_i > 0 else 0 for p_i in p)
    
    # Step 3: Calculate the summation term
    summation_term = 0
    for i in range(nb_classes):
        for x_i in range(nb_seq + 1):
            binom_coeff = math.comb(nb_seq, x_i)
            term = binom_coeff * (p[i] ** x_i) * ((1 - p[i]) ** (nb_seq - x_i))
            log_x_i_fact = math.log(math.factorial(x_i)) if x_i > 0 else 0
            summation_term += term * log_x_i_fact
    
    # Combine all terms
    entropy = -log_n_fact + entropy_term + summation_term
    
    return entropy

def calculate_exp_entrop_sat(alignment:object)->float:
    """ This fucntion calculate the expected entropy under saturation. 
    this is done by looking at the frequency of the four bases in the alignment
    It is expected that all site of the alignment are neutrally evolving

    :param alignemnt:The alignment used to calculate the expeted entropy
    :return: the expected entropy. 
    """
    nb_seq = len(alignment)
    
    # Calculate overall nucleotide frequencies
    overall_probs, overall_counts = get_over_all_nuc_freq(alignment)

    # Calculate expected entropy under full saturation
    H_X_expected = calculate_entropy(nb_seq, list(overall_probs.values()))
 
    return H_X_expected
    
def calculate_entropy_of_alignment(alignment:object)->tuple:
    """This functiona calculate the observed entropy 
    for each column of the alignemnt.

    :param alignment: multiple sequence alignment
    :return (average_entropy, site_entropies): tuble returning the average entropy 
                                               calculated form each site and the 
                                               list of entropy per sites.
    """
    num_sites = len(alignment[0])
    nb_seq = len(alignment)
    site_entropies = []
    for site_idx in range(num_sites):
        site_counts = {"A":0, "T":0, "G":0, "C":0}
        for seq in alignment:
            nucleotide = seq[site_idx]
            if nucleotide in site_counts:
                site_counts[nucleotide] += 1
            else:
                site_counts[nucleotide] = 1
        total_nucleotides = sum(site_counts.values())
        site_probs = {nuc: count / total_nucleotides for nuc, count in site_counts.items()}
        site_entropy = calculate_entropy(nb_seq, list(site_probs.values()))

       # site_entropy = calculate_information_content(site_counts, overall_probs)
        site_entropies.append(site_entropy)
    
    average_entropy = sum(site_entropies) / num_sites
    return average_entropy, site_entropies

def get_over_all_nuc_freq(alignment:object)->tuple:
    """ This function calculate the expected frequency of each nucleotude 
    class in each site. This is done by  calculating the frequency form the 
    whole alignement . Only parsimany informative sites are used for this 
    to remove bias caused by different in evolutionaryrate or difference in 
    number of sites whihc coud affect result.

    :param alignment: the alignemnt used to calculate the expected frequencies
    :return (overall_probs, overall_counts): two dictionary storing count and 
                                             frequency for each nuc classe
    """
    # Calculate overall nucleotide frequencies
    nb_sites = len(alignment[0])
    nb_seq = len(alignment)
    i = 0
    overall_nuc_class = {'A':0, 'C':0,'T':0,'G':0}
    while i < nb_sites:
        column = alignment[:,i:i+1]
        nuc_class = {'A':0, 'C':0,'T':0,'G':0}
        for nuc in column:
            if nuc.seq == "-" or nuc.seq == "N":
                break
            nuc_class[nuc.seq] += 1
        tag_inf = True
        for nuc, val in nuc_class.items():
            if val == nb_seq:
                tag_inf = False # in this the nuc class is the on one represented 
                                   # in the column and so this position is not informative
                break
        if tag_inf:
            for nuc , val in nuc_class.items():
                overall_nuc_class[nuc] += val
        i += 1
    print(overall_nuc_class)
    total_nucleotides = sum(overall_nuc_class.values())
    overall_probs = {nuc: count / total_nucleotides if total_nucleotides > 0 else 0 for nuc, count in overall_nuc_class.items()}
    return overall_probs, overall_nuc_class

def main(align_file:str, out_file:str)-> None:
    """ The main function of the script
    
    :param align_file: path to the alignment file
    :param out_file: file where saturation infroantio are saved
    """
    # load alignemnt
    alignment = load_align(align_file)

    # Calculate expeted entropy of the alignment
    sat_expected_entropy = calculate_exp_entrop_sat(alignment)
    print(f"Expected Entropy under saturation H(X): {sat_expected_entropy}")

    # Calculate observed entropy of the alignment
    mean_observed_entropy, site_entropies = calculate_entropy_of_alignment(alignment)
    print(f"Mean Observed Entropy: {mean_observed_entropy}")

    # Using the Stats library, compute t-statistic and p-value
    t_stat, p_val = stats.ttest_1samp(a=site_entropies, popmean = sat_expected_entropy)
    print("t-statistic = " + str(t_stat))  
    print("p-value = " + str(p_val)) 
    
    if p_val < 0.05:
        print("The observed entropy is significantly different from the expected entropy under full saturation.")
    else:
        print("The observed entropy is not significantly different from the expected entropy under full saturation.")

    # save output
    with open(out_file,"w") as file_handler:
        file_handler.write(f"exp_entrop\tobs_entrop\tt_stat\tpval\n")
        file_handler.write(f"{sat_expected_entropy}\t{mean_observed_entropy}\t{t_stat}\t{p_val}\n")



########################################################################################
########### Main script
########################################################################################


parser = argparse.ArgumentParser(description='This script test whether the overall sites within an alignemnt show evidence for saturation or not')
parser.add_argument('--mult',type=str, help='multiple alignment in fasta format')
parser.add_argument('--out', type=str, help='path to the output file')

args = parser.parse_args()
main(args.mult, args.out)
