#!/user/bin/env python

import pandas as pd
import argparse
import csv
from argparse import RawDescriptionHelpFormatter


# Author; Felicia New, 2018
# Filter RPKM output file based on the coverage column

def getOptions():
    """ Function to pull in command line arguments """
    description="""This script is used to filter the RPKM output based on the coverage across the gene"""

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    
    parser.add_argument("-i","--input_file",dest="input",action="store",required=True,help="Input RPKM file to be filtered [Required].", metavar="INPUT")
    parser.add_argument("-t","--threshold",dest="threshold",action="store",type=int,required=True,default=80,help="Filtering threshold; default: 80 for >80%, [Required]",metavar="THRESHOLD")
    parser.add_argument("-o","--out", dest="out",action="store",required=True,help="Output file for filtered counts [Required]", metavar="OUT")
    args=parser.parse_args()
    return(args)


def filterCounts(args):
    """ Function to read in RPKM file and filter based on the third column and the given threshold cut off """
    cutoff = args.threshold
    df = pd.read_csv(args.input, sep=",")
    filt = df.loc[df[df.columns[2]] > cutoff]
    return(filt)

def writeOutput(args, filt):
    final = filt.drop(filt.columns[2], axis=1)
    final.to_csv(args.out, sep=",", index=False)

def main():
   args = getOptions()
   filt = filterCounts(args)
   writeOutput(args, filt)

if __name__=='__main__':
    main()

