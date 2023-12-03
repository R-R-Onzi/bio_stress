import pandas as pd
from os import walk
import argparse

def main(args, fantom_files: list):

    results_df = pd.read_csv(args.results_file, delimiter=",",header=True)
    
    paths_df = pd.read_csv(args.path_file, delimiter=",",header=True)

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

    parser.add_argument('--path_file', metavar='path_file', required=True, help='path file')

    parser.add_argument('--phantom_folder', metavar='phantom_folder', required=True, help='phantom folder')
    
    parser.add_argument('--results_file', metavar='results_file', required=True, help='results file')
    
    args = parser.parse_args()
    
    fantom_files = []
    
    for (dirpath, dirnames, filenames) in walk(args.phantom_folder):
        for filename in filenames:
            list_rez: list = []

            if ".csv" not in filename or "p1" not in filename or "," in filename:
                continue

            filename_split = filename.split(sep="@")
                
            list_rez.append(filename_split[0] + "@")
            list_rez.append(filename_split[1].split(".")[0])
            list_rez.append("." + filename_split[1].split(".")[1])
            fantom_files.extend([list_rez])
        
    main(args, fantom_files)
