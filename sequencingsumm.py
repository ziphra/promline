import pandas as pd
import os
import subprocess
import glob
import sys


fast5 = sys.argv[1]
print(fast5)
parent_folder = os.path.dirname(fast5)
pattern = os.path.join(parent_folder, 'sequencing_summary*')
summ = glob.glob(pattern)
print(summ)


if summ:
    target_file = summ[0]
    print(target_file) # Output: /path/to/sequencing_summary.txt (or similar)

    # reading files
    summary = pd.read_csv(summ[0],sep='\t')
    qscore = pd.read_csv('h_qscore.txt',sep='\t')
    rl = pd.read_csv('h_read.txt',sep='\t')
    
    # calculate read length from read sequence
    rl['sequence_length_template']  = rl['read'].str.len()

    # using merge function
    output1 = pd.merge(summary, qscore, on='read_id')

    # remove read sequence
    rl.drop("read",axis=1,inplace=True)

    # merge again 
    output2 = pd.merge(output1, rl, on='read_id')

    # save file
    output2.to_csv("presequencing_summary.txt",index=False,sep='\t')

else:
    print("No summary file = no QC")
