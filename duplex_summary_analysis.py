#!/usr/bin/python
# Heather Borgard (hborgard@ucsc.edu)
# duplex_summary_analysis.py

import os, sys, time
import numpy as np
import gzip
import pandas as pd
import argparse

########################################################################
# Main
# Here is the main program
########################################################################

def main(myCommandLine=None):

    t0 = time.time()

    #Parse the inputs arguments allowing multiple files and optional -n sample name 
    parser = argparse.ArgumentParser()
    parser.add_argument('nfile', nargs='+')
    parser.add_argument("-n", "--name", help='the sample name')

    #Parse the arguments
    args = parser.parse_args()

    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    #option to unzip the file if needed 
    for inFile in args.nfile:
        if inFile.endswith('.gz'):
            file = gzip.open(inFile,'rt')
        else:
            file = open(inFile, 'r')
    
        #convert txt file to csv with whitespace deliminter to use dataframes 
        df=pd.read_csv(file, delim_whitespace=True )
        
        #print(list(df.columns))
        
        #drop unnecessary columns for easier processing 
        #df=df.drop(columns=['filename', 'read_id', 'run_id', 'channel', 'mux', 'start_time', 'duration', 'template_start', 'template_duration'])
        
        #create bins on qscore
        q10=df[df['mean_qscore_template'] >= 10]
        q20=df[df['mean_qscore_template'] >= 20]
        q30=df[df['mean_qscore_template'] >= 30]
        q40=df[df['mean_qscore_template'] >= 40]
        q50=df[df['mean_qscore_template'] >= 50]
        #print(quals10.value_counts())
        
        #determine bases (using length template for quality bins)
        bases=df['sequence_length_template']
        bases10=q10['sequence_length_template']
        bases20=q20['sequence_length_template']
        bases30=q30['sequence_length_template']
        bases40=q40['sequence_length_template']
        bases50=q50['sequence_length_template']

        #sum up total bases for each quality bin
        total_bases_10 = np.sum(bases10)
        total_bases_20 = np.sum(bases20)
        total_bases_30 = np.sum(bases30)
        total_bases_40 = np.sum(bases40)
        total_bases_50 = np.sum(bases50)

        #get other stats that arent quality binned 
        read_length = sorted(bases, reverse=True)
        total_bases = np.sum(bases)
        total_gigabases = round(total_bases / 1E9, 2)
        target = total_bases / 2.0

        
        # these stats are based on the human genome (3.3E9 bases per genome)
        coverage = round(total_bases / 3.3E9, 2)
        lt100 = round(sum([i for i in read_length if i >= 100000]) / 3.3E9, 2)
        lt200 = round(sum([i for i in read_length if i >= 200000]) / 3.3E9, 2)
        lt300 = round(sum([i for i in read_length if i >= 300000]) / 3.3E9, 2)
        lt400 = round(sum([i for i in read_length if i >= 400000]) / 3.3E9, 2)
        lt500 = round(sum([i for i in read_length if i >= 500000]) / 3.3E9, 2)
        lt1000 = round(sum([i for i in read_length if i >= 1000000]) / 3.3E9, 2)

        #get coverage for different qscores
        q10c = round(total_bases_10 / 3.3E9, 2)
        q20c = round(total_bases_20 / 3.3E9, 2)
        q30c = round(total_bases_30 / 3.3E9, 2)
        q40c = round(total_bases_40 / 3.3E9, 2)
        q50c = round(total_bases_50 / 3.3E9, 2)

    file.close()

    #get filename based on -n optional argument or naming convention
    if args.name:
        filename = args.name
    else:
        filename = sys.argv[1:]
    
    # get number of 1 Mb+ reads
    num1000 = len([i for i in read_length if i >= 1000000])

    #print the headers of the output 
    print('Sample\tread_N50\tGb\tcoverage\t100kb+\tq10\tq20\tq30\tq40\tq50')


    n50 = 0
    counter = 0
    for item in read_length:
        n50 += item
        counter += 1
        if n50 >= target:
            print(filename, '\t', read_length[counter], '\t', total_gigabases, '\t', coverage, '\t', lt100, '\t', q10c, '\t', q20c, '\t', q30c, '\t', q40c, '\t', q50c, file=sys.stdout)
            break

    print('\ntotal time for the program %.3f' % (time.time()-t0), file=sys.stderr)

if (__name__ == '__main__'):
    main()
    raise SystemExit

