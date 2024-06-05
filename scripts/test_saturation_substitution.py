import math
import numpy as np
from scipy.stats import ttest_ind
from Bio import AlignIO
import argparse 

def load_align(mult_fasta:str)->object:
    """
    Load multiple peptide sequence alignment

    :param mult_fasta:  multiple alignement file in fasta format
    :return: Alignment object
    """
    alignment = None
    with open(mult_fasta) as file_handler:
        alignment = AlignIO.read(file_handler, "fasta")
    return alignment



def calculate_H_X(n, p):
    """
    """
    k = len(p)
    
    # Step 1: Calculate -log(n!)
    log_n_fact = math.log(math.factorial(n))
    
    # Step 2: Calculate -n * sum(p_i * log(p_i))
    entropy_term = -n * sum(p_i * math.log(p_i) if p_i > 0 else 0 for p_i in p)
    
    # Step 3: Calculate the summation term
    summation_term = 0
    for i in range(k):
        for x_i in range(n + 1):
            binom_coeff = math.comb(n, x_i)
            term = binom_coeff * (p[i] ** x_i) * ((1 - p[i]) ** (n - x_i))
            log_x_i_fact = math.log(math.factorial(x_i)) if x_i > 0 else 0
            summation_term += term * log_x_i_fact
    
    # Combine all terms
    H_X = -log_n_fact + entropy_term + summation_term
    
    return H_X

def calculate_information_content(counts, overall_probs):
    """
    """
    total_counts = sum(counts.values())
    
    # Calculate the probability of this site under the multinomial distribution
    multinomial_prob = math.factorial(total_counts)
    for nucleotide, count in counts.items():
        multinomial_prob *= overall_probs[nucleotide] ** count / math.factorial(count)
    
    # Calculate the information content
    information_content = -math.log2(multinomial_prob)
    return information_content

def calculate_entropy_of_alignment(alignment, overall_probs):
    """
    """
    num_sites = len(alignment[0])
    
    site_entropies = []
    for site_idx in range(num_sites):
        site_counts = {}
        for seq in alignment:
            nucleotide = seq[site_idx]
            if nucleotide in site_counts:
                site_counts[nucleotide] += 1
            else:
                site_counts[nucleotide] = 1
        site_entropy = calculate_information_content(site_counts, overall_probs)
        site_entropies.append(site_entropy)
    
    average_entropy = sum(site_entropies) / num_sites
    return average_entropy, site_entropies

def simulate_entropies(n, overall_probs, num_sites, num_simulations=1000):
    """
    """
    simulated_entropies = []
    for _ in range(num_simulations):
        simulated_alignment = []
        for _ in range(n):
            simulated_seq = ''.join(np.random.choice(list(overall_probs.keys()), p=list(overall_probs.values()), size=num_sites))
            simulated_alignment.append(simulated_seq)
        
        _, simulated_site_entropies = calculate_entropy_of_alignment(simulated_alignment, overall_probs)
        simulated_entropies.extend(simulated_site_entropies)
    
    return simulated_entropies

def test_saturation(alignment):
    """
    """
    n = len(alignment)
    num_sites = len(alignment[0])
    
    # Calculate overall nucleotide frequencies
    overall_counts = {}
    for seq in alignment:
        for nucleotide in seq:
            if nucleotide in overall_counts:
                overall_counts[nucleotide] += 1
            else:
                overall_counts[nucleotide] = 1
    total_nucleotides = sum(overall_counts.values())
    overall_probs = {nuc: count / total_nucleotides for nuc, count in overall_counts.items()}
    
    # Calculate expected entropy under full saturation
    H_X_expected = calculate_H_X(n, list(overall_probs.values()))
    print(f"Expected Entropy H(X): {H_X_expected}")
    
    # Calculate observed entropy of the alignment
    observed_entropy, site_entropies = calculate_entropy_of_alignment(alignment, overall_probs)
    print(f"Observed Entropy: {observed_entropy}")
    
    # Simulate entropies to create a distribution
    simulated_site_entropies = simulate_entropies(n, overall_probs, num_sites)
    
    # Perform a t-test
    t_stat, p_value = ttest_ind(site_entropies, simulated_site_entropies)
    print(f"t-statistic: {t_stat}, p-value: {p_value}")
    
    if p_value < 0.05:
        print("The observed entropy is significantly different from the expected entropy under full saturation.")
    else:
        print("The observed entropy is not significantly different from the expected entropy under full saturation.")
    
    return H_X_expected, observed_entropy, t_stat, p_value


def main(align_file, out_file)-> None:
    
    alignment = load_align(align_file)
    H_X_expected, observed_entropy, t_stat, p_value = test_saturation(alignment)

    with open(out_file,"w") as file_handler:
        file_handler.write(f"exp_entrop\tobs_entrop\tt_stat\tpval\n")
        file_handler.write(f"{H_X_expected}\t{observed_entropy}\t{t_stat}\t{p_value}\n")




########################################################################################
########### Main script
########################################################################################


parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--mult',type=str, help='multiple alignment in fasta format')
parser.add_argument('--out', type=str, help='path to the output file')
args = parser.parse_args()

main(args.mult, args.out)
