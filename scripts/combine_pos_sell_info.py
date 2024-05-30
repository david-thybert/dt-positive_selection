import argparse 
import json
import pandas as pd
import statsmodels.stats.multitest as ssm

#from statsmodels import stats

def fetch_pos_sel_info(json_file:str)->dict:
    """
    Retrieve the positive selction information fomr the json file

    :param json_file: JSON file from positive selction analysis
    :return: dictionary of results
    """

    with open(json_file) as json_handler:
        file_contents = json_handler.read()
        
    parsed_json = json.loads(file_contents)

    # gee if is contained in the file name
    gene_id = json_file.split("/")[-1].split(".nuc")[0]
    result = {}
    for (branch,values) in parsed_json["branch attributes"]['0'].items():
        if not "original name" in values:
            continue
        species_name = values["original name"].split("|")[-1]
        pval = values['Uncorrected P-value']
        pval_corr = values['Corrected P-value']
        rate_class = values['Rate classes']
        lrt = values['LRT']
        omega_ratio_base = values['Baseline MG94xREV omega ratio']
        base = values['Baseline MG94xREV']
        if not pval is None:
            result[species_name] = [gene_id, lrt, pval, pval_corr, rate_class, omega_ratio_base, base]
    return result

def _create_data_frame(pos_sel_branches:dict)->dict:
    """
    This function create a dataframe fromt the matrix of postive selection informaton associated
    to each branch/species in the input dictionary.

    :param pos_sel_branches: dictionary storing the information matrix for each species/branch tested
    :return: a dictionary with data frame instead of matrix associated to each branch/species tested.
    """
    result = {}
    for branche, pos_val in pos_sel_branches.items():
        df = pd.DataFrame(pos_val,  columns =  ["gene_id", "lrt", "pval", "pval_corr", 
                                                "rate_class", "omega_ratio_base", "base"])
        result[branche] = df
    return result

def multitetesting_correction(pos_sel_branches_df:dict, method:str="fdr_bh")->dict:
    """
    This function performed a multitesting correction
    
    :param pos_sel_branches_df: the dictionary that store for each branch/species a 
                                dataframe with the positicve selection results 
    :param method: the method to be used for multitesting correction
    :return: dictionary of that store the postive selction data for each branch/species
             including the adjusted pvalues.
    """
    for branch, values_df in pos_sel_branches_df.items():
        pval = values_df["pval_corr"]
        rej, pval_adj, alphasidak, alphacBonf = ssm.multipletests(pval, method=method)
        values_df["pval_adj"] = pval_adj
    return pos_sel_branches_df

def main(files:str, pref_out:str)->None:
    """
    The main fucntion of the script

    :param files: string that list all json files to be parsed separated by a white space
    :param pref_out: prefix used for the file name
    """
    lst_jsons = files.split()
    pos_sel_branches = {}
    for json in lst_jsons:
        pos_sel = fetch_pos_sel_info(json)
        for (species, val) in pos_sel.items():
            if not species in pos_sel_branches:
                pos_sel_branches[species] = []
            pos_sel_branches[species].append(val)
    
    #converting array to pandas df
    pos_sel_branches_df = _create_data_frame(pos_sel_branches)
    
    # mutli testtin correction
    pos_sel_branches_mlt = multitetesting_correction(pos_sel_branches_df)
   
    #save the file per branche
    for branche, df in pos_sel_branches_mlt.items():
        file_name = f"{pref_out}_{branche}.possel"
        df.to_csv(file_name,sep="\t")


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--files',type=str, help='list of json file to integrate')
parser.add_argument('--prefix_out', type=str, help='prefix for outfiles')
args = parser.parse_args()

main(args.files, args.prefix_out)

