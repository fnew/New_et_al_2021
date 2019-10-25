#!/usr/bin/env python

from __future__ import division
import argparse
import logging
import os
import csv
import collections
import numpy as np
from argparse import RawDescriptionHelpFormatter
csv.field_size_limit(1000000000)

# Author: Felicia New, 2017


def getOptions():
    """ Function to pull in arguments """

    description="""This script can be used to calculates coverage (RPKM) using a BAM file and the output from SAMTools depth.\n
    RPKM = Reads per kilobase per million, (Mortazavi, 2008)"""



    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-f1", "--fastq1", dest="fastq1", action='store',type=int, required=True, help="Length of fastq file 1 (wc -l) [Required]", metavar="FASTQ1")     
    parser.add_argument("-f2", "--fastq2", dest="fastq2", action='store', type=int, required=False, help="Length of fastq file 2 (wc -l) [Optional -- if PE sequencing]", metavar="FASTQ2")     
    parser.add_argument("-s", "--stat", dest="sname", action='store', required=True, help="Output from samtools idxstats [Required]", metavar="STAT_FILE")
    parser.add_argument("-d", "--depth", dest="depth", action='store', required=True, help="Output from samtools depth [Required]", metavar="DEPTH")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Out File", metavar="OUT_FILE")
    args = parser.parse_args()
    return(args)


def setLogger(fname, loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')


# SAM functions
#def read_bam(args):
#    """ Function to read BAM file to get number of mapped reads """
#    logging.info("Reading BAM file '%s'." % args.bname)
#    samfile = pysam.AlignmentFile(args.bname, "rb")
#    total_mapped_reads=int(samfile.mapped)
#    return(total_mapped_reads)

# Get total reads per sample
def count_reads(args):
    """ Taking the fastq line counts and dividing by 4 to get the number of reads per sample. This function will sum the R1 and R2 files. """
    logging.info("Getting number of reads per sample")
    len1 = int(args.fastq1) / 4
    len2 = int(args.fastq2) / 4
    total_reads_in_sample = len1 + len2
    return(total_reads_in_sample)


# IDXStats Functions
def read_stats(args):
    """ Reading the IDXStats output file to get region_length and number of reads mapped to each region """
    logging.info("Creating and reading the idxstats file")
    cdict = collections.defaultdict(dict)
    with open(args.sname,'r') as STAT:
        reader = csv.reader(STAT,delimiter='\t')
	for row in reader:
	    name = row[0]
	    length = int(row[1])
	    num_mapped = int(row[2])
	    num_unmapped = int(row[3])
	    start = 0
	    end = length
	    cdict[name]['name'] = name #Name of the gene/reference
            cdict[name]['length'] = length #Length of the gene
	    cdict[name]['num_mapped'] = num_mapped #Number of reads mapped to the gene
	    cdict[name]['num_unmapped'] = num_unmapped
	    cdict[name]['start'] = start
	    cdict[name]['end'] = end
	    cdict[name]['pos'] = {}
        return(cdict)

# MPILEUP functions
def read_mpileup(args,cdict):
    """ Read mpileup and store depth and length into dictionary """
    logging.info("Reading the mpileup file '%s'." % args.depth)
    with open(args.depth, 'r') as DEPTH:
        reader = csv.reader(DEPTH, delimiter='\t', quoting=csv.QUOTE_NONE)
        
	for row in reader:
	    gene = row[0]
	    mpos = int(row[1]) - 1 #mpileups are 1-based
	    mdepth = int(row[2])
            cdict[gene]['pos'][mpos]=mdepth
    for gene in cdict:
	bpcov=0 #create flag for counting bases that have coverage
	for pos in cdict[gene]['pos']:
            bpcov += np.count_nonzero(cdict[gene]['pos'][pos])
	cdict[gene]['bpcov']=bpcov
    return(cdict)

# Coverage Functions
    """ Calculate different coverage metrics: Estimate number of reads in regions, Average per nucleotide coverage (APN),
    Reads per kilobase per million mapped reads (RPKM) and coverage accross the gene: what proportion of bp have reads mapped to them. """
def calc_coverage(args,cdict,total_reads_in_sample):
    logging.info("Calculating Coverage Counts")
    for name in cdict:
	depth = (cdict[name]['num_mapped'])
	
	if depth != 0:
	    adj_length = cdict[name]['length'] / 1000
	    adj_total_reads = total_reads_in_sample / 1000000
	    coverage = (cdict[name]['bpcov'] / cdict[name]['length']) * 100 #Percentage of the gene that has coverage
	    cdict[name]['coverage'] = coverage
	    cdict[name]['rpkm'] = (cdict[name]['num_mapped']) / (adj_length * adj_total_reads) # Calculates reads per kilobase per million reads RPKM from Moretzavi et al.
	    
	else:
	    # if there is no coverage set everything to 0
	    cdict[name]['coverage'] = 0
	    cdict[name]['rpkm'] = 0
    return(cdict)


#Output Functions
def writeOutput(args,cdict):
    """ Brute force methods... """
    logging.info("Writing Output")
    fname=os.path.basename(args.depth)
    header = ['gene_id','rpkm_%s' % fname, 'coverage_%s' % fname]
    with open(args.out, 'wb') as OUT:
	OUT.write(','.join(header) + '\n')
	for name in cdict:
	    OUT.write(','.join(str(name) for name in [name,cdict[name]['rpkm'],cdict[name]['coverage']]) + '\n')



def main():
    """ MAIN Function to execute everything """
    args = getOptions()
    if args.log:                               # Turn on logging if option -g was given
	setLogger(args.log,logging.INFO)

    total_reads_in_sample = count_reads(args)  #Use fastq file to count the number of total reads in a sample
    cdict = read_stats(args) #Read through stats file and create dictionary to sort all info
    read_mpileup(args,cdict)  #Read the mpileup/depth file to get pos coverage info
    calc_coverage(args,cdict,total_reads_in_sample)    # Use info stored in bdict to calculate coverage (APN, RPKM), and other measures
    writeOutput(args,cdict) #Write output to CSV file

if __name__=='__main__':
    main()
    logging.info("Script complete")


