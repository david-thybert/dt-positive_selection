import argparse 
from ete3 import Tree


"""
This script launch codeml and check that there is not errors

command:

    python tag_tree_nodes.py --tree <nwk file> --spe_tag <file listing leafs to tag> --tag_name <"tag_name"> --ancestral <True/False> --out <tagged nwk file>

"""

def load_species_to_tag(species_tag:str)->list:
    """
    Load the list fo species to tag in the phyloegentic tree

    :param species_tag: path the to the file that contain the list of species
    :return: list of species to be tagged
    """
    results = []
    with open(species_tag) as file_handler:
        for line in file_handler:
            results.append(line.strip())
    return results

def tag_only_ancestors(lca:Tree, tag:str)->Tree:
    """ Tag the ancestral node only

    :param lca: last common ancestor node
    :param tag: tag to be added to the node
    return : lca node with the tag
    """
    if tag not in lca.name:
        new_name = lca.name + tag
        lca.add_feature("name", tag)
    return lca

def tag_ancestors(node:Tree, lca:Tree, tag:str)->Tree:
    """
    Tag all the path from the node to the last common ancestor.
    Return the node tagged

    :param node: node starting point to anotated with the tag and all ancestors until lca
    :param lca: last common ancestor node
    :param tag: tag to be added to the node
    :return: node and their ancestro tagged
    """
    curr_node = node

    while curr_node != lca:
        if tag in curr_node.name:# if already tagged, no need to add the tag
            curr_node = curr_node.up
            continue 
        new_name = curr_node.name + tag
        curr_node.add_feature("name", new_name)
        curr_node = curr_node.up

    if curr_node == lca and  tag not in curr_node.name:
        new_name = curr_node.name + tag
        curr_node.add_feature("name", tag)
    return node

def tag_tree(tree:Tree, species_to_tag:list, tag_name:str, anc:bool)->Tree:
    """
    This function tag all leafs deined in the list to tag. If the ancestor 
    tag is true then it will tag all ancestor from each leaf until the LCA 

    :param tree: the tree containing the tags.
    :param species_to_tag: list of species to tag
    :param tag_name: the string corresponding to the tag"
    :param anc: True tag the ancestor and LCA , False tag just the leafs.
    """
    leaf_nodes = []
    species_with_node = []
    for species in species_to_tag:
        print(species)
        for leaf in tree:
            if species == leaf.name.split("|")[-1]:
                print(leaf.name)
                leaf_nodes.append(leaf)
                species_with_node.append(species)
                break
    print(leaf_nodes[1:])
    if anc == 1 and len(species_to_tag) > 1:# if we want to tag all ancestors until lca we should have more than 1 species
        lca_node = leaf_nodes[0].get_common_ancestor(leaf_nodes[1:])
        for node in leaf_nodes:
            tag_ancestors(node, lca_node, tag_name)
    elif anc == 2 and len(species_to_tag) == 2:
        lca_node = leaf_nodes[0].get_common_ancestor(leaf_nodes[1:])
        tag_only_ancestors(lca_node, tag_name)
    elif anc == 0:
        for node in leaf_nodes:
            new_name = node.name + tag_name
            node.add_feature("name", new_name)
    else:
        raise Exception(f"inconcistancy between parameters anc({anc}) and number of species({len(species_with_node)})\n if anc == 1 thene the numebr of species need to be at least 2 \n if anc == 2 then the number of species need to be only 2")
    return tree

def main(tree_file:str, species_tag:str, tag_name:str, anc:int, out_file:str)->None:
    """
    Main fuunction of the script

    :param tree_file: path to the tree file in nwk format
    :param species_tag: path to the file with the name of species to tag. 
                        One species per line.
    :param tag_name: the tag name
    :param anc: True tag ancestors until the LCA between all species to tag.
    :param out_file: path to the outfile where to store the tagged tree in nwk format.
    """
    tree = Tree(tree_file, format = 1)
    species_to_tag = load_species_to_tag(species_tag)
    tagged_tree = tag_tree(tree, species_to_tag, tag_name, anc)
    tagged_tree.write(outfile=out_file, format=1)


########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='This script tag a branch of nwk tree to be used for phylogenetic analyses')
parser.add_argument('--tree',type=str, help='nwk file describing the phylogeny')
parser.add_argument('--spe_tag', type=str, help='file that contains the list of species to tag in the tree')
parser.add_argument('--out', type=str, help='path to the tagged tree')
parser.add_argument('--tag_name', type=str, help='name of the tag for the forgrand analysis')
parser.add_argument('--ancestral', type=int, default=False, help='tag defining if and how cancestral nod will be taged (0: no ancestor tagged, 1: leafs and ancestor, 2: only ancestor. the tag 2 need to ahve only two species defined in the config file)')

args = parser.parse_args()
main(args.tree, args.spe_tag, args.tag_name, args.ancestral, args.out)
