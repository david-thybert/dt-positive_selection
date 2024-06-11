import argparse 
import json
import pandas as pd
import statsmodels.stats.multitest as ssm

#from statsmodels import stats

def map_inputs(tsvs:list, sat_substs:list)-> dict:
    """
    This function map together json file and saturation of substitution file

    :param jsons: list of json fles
    :param sat_substs: list of saturaiton of substitution files
    :return: Dictionary that map json file tot he corresponding saturations substtitution file
    """
    dico_result = {}
    for tsv in tsvs:
        id_tsv = tsv.split("/")[-1].split(".pml")[0]
        for sta_subst in sat_substs:
            id_sta_subst = sta_subst.split("/")[-1].split(".nuc")[0]
            if id_tsv == id_sta_subst:
                dico_result[id_tsv] = [tsv, sta_subst]
    return dico_result

def fetch_pos_sel_info(gene_id:str, tsv_file:str, sat_subst:str)->dict:
    """
    Retrieve the positive selction information fomr the json file

    :param json_file: JSON file from positive selction analysis
    :return: dictionary of results
    """
    lst_val = []
    with open(sat_subst) as sat_subst_handler:
        for line in sat_subst_handler:
            if "exp_entrop" in line:
                continue
            lst_val = line.split()
            break

    with open(tsv_file) as file_handler:
        file_contents = file_handler.readlines()
    
    parsed_tsv = file_contents[1].split()
    values = []
    for val in parsed_tsv[1:]:
        values.append(float(val))
    result =  parsed_tsv[0:1] + values + [float(lst_val[-1]), float(lst_val[1]), float(lst_val[0])]
    return result

def _create_data_frame(pos_sel_branches:list)->dict:
    """
    This function create a dataframe fromt the matrix of postive selection informaton associated
    to each branch/species in the input dictionary.

    :param pos_sel_branches: dictionary storing the information matrix for each species/branch tested
    :return: a dictionary with data frame instead of matrix associated to each branch/species tested.
    """

    df = pd.DataFrame(pos_sel_branches,  columns =  ["gene_id", "lrt", "pval", "alt_w_0", 
                                                    "alt_w_1", "alt_w_2a", "alt_w_2b", "null_w_0", 
                                                    "null_w_1", "null_w_2a", "null_w_2b", "Pval_no_sat", 
                                                    "obs_entropy", "exp_entropy"])
    return df

def multitetesting_correction(pos_sel_df:object, method:str="fdr_bh")->dict:
    """
    This function performed a multitesting correction
    
    :param pos_sel_branches_df: the dictionary that store for each branch/species a 
                                dataframe with the positicve selection results 
    :param method: the method to be used for multitesting correction
    :return: dictionary of that store the postive selction data for each branch/species
             including the adjusted pvalues.
    """
    pval = pos_sel_df["pval"]
    rej, pval_adj, alphasidak, alphacBonf = ssm.multipletests(pval, method=method)
    pos_sel_df["pval_adj"] = pval_adj
    return pos_sel_df

def main(files:str, file_sat_subst:str, pref_out:str)->None:
    """
    The main fucntion of the script

    :param files: string that list all json files to be parsed separated by a white space
    :param file_sat_subst: the list fo file with substitution saturaiton info
    :param pref_out: prefix used for the file name
    """
    tsvs = files.split()
    subst_sats = file_sat_subst.split()
 
    dico_input = map_inputs(tsvs, subst_sats)

    pos_sel_vals = []
    for id, files in dico_input.items():
        tsv = files[0]
        sat_subst = files[1]
        val = fetch_pos_sel_info(id, tsv, sat_subst)
        pos_sel_vals.append(val)

    #converting array to pandas df
    print(pos_sel_vals)
    pos_sel_val_df = _create_data_frame(pos_sel_vals)
    print(pos_sel_val_df)
    # mutli testtin correction
    pos_sel_val_mlt = multitetesting_correction(pos_sel_val_df)
   
    #save the file per branche
    file_name = f"{pref_out}.possel"
    pos_sel_val_mlt.to_csv(file_name, sep="\t")


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--files_possel',type=str, help='list of json file to integrate')
parser.add_argument('--files_sat_subst',type=str, help='list of files for saturation_substitution')
parser.add_argument('--prefix_out', type=str, help='prefix for outfiles')
args = parser.parse_args()

main(args.files_possel, args.files_sat_subst, args.prefix_out)

