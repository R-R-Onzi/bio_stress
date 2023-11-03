import pandas as pd
import argparse


def main(file: str):

    df = pd.DataFrame(pd.read_excel(file)) 

    df.to_csv(f"{file.split('.')[0]}.csv")

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', metavar='file path', required=True,
    help='file folder')
    
    args = parser.parse_args()
    
    main(args.file)
