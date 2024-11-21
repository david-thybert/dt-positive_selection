import sys
import argparse
import pandas as pd



def is_sequential_pos(row:object,number:int,windows:int)-> bool:
    """
    """
    pos_string = row["position"]
    positions = pos_string.split("|")
    last_pos = -1
    nb = 0
    i=0
    while i < len(positions):
        p = positions[i]
        pos = int(p.split(":")[0])
        lst_pos = [pos]
        j = i + 1
        while j < len(positions):
            p_next = positions[j]
            pos_next = int(p_next.split(":")[0])
            if pos_next - pos <=windows:
                lst_pos.append(pos_next)
            else:
                if len(lst_pos)>=number:
                    return True
                break
            j = j +1
        i = i +1
    return False


def main(possel:str, threshold:float, nb:int, window:int, out:str)-> None:
    """
    """
    df = pd.read_csv(possel, index_col=0, sep="\t")
    print(df.head())
    df['pval_adj'] = df['pval_adj'].astype(float)
    significant_df = df[df["pval_adj"] <= threshold]

    filtered_df = pd.DataFrame()
    for index, row in significant_df.iterrows():
        if not is_sequential_pos(row, nb, window):
            filtered_df.append(row, ignore_index=True)
    filtered_df.to_csv(out, sep="\t", index=False)





###########################################################################
#########  Main programm
###########################################################################

parser = argparse.ArgumentParser(description='Script formating the data to be handle by the positvie selction pipeline')
parser.add_argument('--possel',type=str, help='possel file')
parser.add_argument('--thr',type=float, default=0.05, help='adj pval threshold to filter')
parser.add_argument('--nb',type=int, help='possel file')
parser.add_argument('--wind',type=int, default=0.05, help='adj pval threshold to filter')
parser.add_argument('--out', type=str, help='outfile')

args = parser.parse_args()
main(args.possel, args.thr, args.out)
