
import argparse
import os

def main(dir1:str, dir2:str)->None:
    """
    """
    files_1 = os.listdir(dir1)
    files_1_set = set()
    for file in files_1:
        file.split(".")[0]
        files_1_set.add(file.split(".")[0])
    
    files_2 = os.listdir(dir2)
    files_2_set = set()
    for file in files_2:
        file.split(".")[0]
        files_2_set.add(file.split(".")[0])
    
    for elem in files_1_set:
        if not elem in files_2_set:
            print(elem)

########################################################################################
########### Main script
########################################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--dir1',type=str, help='path to first directory')
parser.add_argument('--dir2', type=str, help='path to second directory')

args = parser.parse_args()
main(args.dir1, args.dir2)


