import pandas as pd
from os import walk
from collections import defaultdict

listis = [] 
result_df = pd.DataFrame()
for (dirpath, dirnames, filenames) in walk("GSE235915_counts_lists"):
    listis.extend(filenames)
    break
i = 0
for files in listis:
    df_trips =  pd.read_csv(f'GSE235915_counts_lists/{files}', delimiter="\t",header=None)
    df_trips = df_trips.iloc[4:]
    if(i==0):
        f_column = df_trips.iloc[:, [0]]
        result_df = pd.concat([result_df,f_column], axis = 1)
    f_column = df_trips.iloc[:, [1]]
    result_df = pd.concat([result_df,f_column], axis = 1)
    i+=1
result_df.to_csv("GSE235915.txt", index=False, header=None)
