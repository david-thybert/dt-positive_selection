import sys
import argparse
import json
from pprint import pprint


def load_conf(conf:str)->dict:
    """
    """
    result = {}
    with open(conf) as f:
        result = json.load(f)
        print(result)
    return result

   

def main(dir_ref:str, conf:str, out:str)->None:
    """
    """
    conf = load_conf(conf)






###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--dir',type=str, help='file describing the orhtology relationship between genes ')
parser.add_argument('--conf', type=int, help='path to the directory contianing the CDS files')
parser.add_argument('--out', type=str, help='path to the outdirectory')

args = parser.parse_args()
main(args.dir, args.conf,args.out)

