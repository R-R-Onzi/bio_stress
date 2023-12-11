import pandas as pd
from os import walk
import argparse
from collections import defaultdict
import networkx as nx
import matplotlib
import matplotlib.pyplot
import csv

def main(args, fantom_files: list):

    results_df = pd.read_csv(args.results_file, delimiter=",",header=None) 
    paths_df = pd.read_csv(args.path_file, delimiter=",",header=None) # go and genes

    results_df = results_df.iloc[1:] #remove first row
    paths_df = paths_df.iloc[2:,1:]
    
    results_df = results_df.iloc[(-results_df[5].astype(float).abs()).argsort()] # sort by biggest nes

    uniques = []

    for j in range(len(paths_df)):
        for a in paths_df.iloc[j, 1].split(","):
            if(a.replace(' ', '') not in uniques):

                uniques.append(a.replace(' ', ''))

    with open(f"{args.results_file.replace('.csv', '').split('/')[-1]}_unique_genes.csv", 'a+') as f:
        write = csv.writer(f)

        write.writerow(uniques)



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
