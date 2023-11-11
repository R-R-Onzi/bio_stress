import pandas as pd
from os import walk
import argparse


def main(folder: str, listis: list):
    listis.sort()
    result_df = pd.DataFrame()

    i = 0
    for file in listis:
        df_trips =  pd.read_csv(f'{folder}/{file}', delimiter="\t",header=None)
        df_trips = df_trips.iloc[4:]
        if(i==0):
            f_column = df_trips.iloc[:, [0]]
            result_df = pd.concat([result_df,f_column], axis = 1)
        f_column = df_trips.iloc[:, [1]]
        result_df = pd.concat([result_df,f_column], axis = 1)
        i+=1

    result_df.to_csv(f"{folder.split(sep="_")[0]}.txt", index=False, header=None)


if __name__  == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', metavar='path', required=True,
    help='folder points')
    
    args = parser.parse_args()
    
    listis = []
    
    for (dirpath, dirnames, filenames) in walk(args.folder):
        listis.extend(filenames)
        break
    
    main(args.folder, listis)
