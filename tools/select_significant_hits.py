import sys
import argparse
import pandas as pd

def main(possel:str, threshold:float, out:str)-> None:
    """
    """
    df = pd.read_csv(possel, index_col=0, sep="\t")
    print(df.head())
    df['pval_adj'] = df['pval_adj'].astype(float)
    filteres_Df = df[df["pval_adj"] <= threshold]
    filteres_Df.to_csv(out, sep="\t", index=False)




###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--possel',type=str, help='possel file')
parser.add_argument('--thr',type=float, default=0.05, help='adj pval threshold to filter')
parser.add_argument('--out', type=str, help='outfile')

args = parser.parse_args()
main(args.possel, args.thr, args.out)
