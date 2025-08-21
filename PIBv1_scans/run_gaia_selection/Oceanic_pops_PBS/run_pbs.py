#!/usr/bin/env python3

import os
import pathlib
import re
import itertools
import allel
import numpy as np
import pandas as pd
from glob import glob
from numpy import asarray
from numpy import savetxt

w_size, w_start, w_stop, w_step = 20, 0, None, 5

pop1_files = sorted(glob("pops_oceania/pop1_*.txt"))

file_name = "pops_oceania/pop2_EAS.txt"
with open(file_name,'r') as f:
    pop2 = f.read().split("\n")[:-1]

file_name = "pops_oceania/pop3_AF-EU.txt"
with open(file_name,'r') as f:
    pop3 = f.read().split("\n")[:-1]

# if True: (and remove the for loop)
# vcf_filename = 'PIBv1_chr6_VQSRpassMissingness0.05AllGTmasks_MAFannotated.SNP_only.vcf.gz'
vcf_files = sorted(glob('../../vcf_files/PIBv1_chr*_final_phased.SNP_only.vcf.gz'))

# for vcf_filename in vcf_files:
def run_vcf_file(vcf_filename, pop1_filename):
    with open(pop1_filename,'r') as f:
        pop1 = f.read().split("\n")[:-1]
    callset_pop1 = allel.read_vcf(vcf_filename, samples=pop1)
    callset_pop2 = allel.read_vcf(vcf_filename, samples=pop2)
    callset_pop3 = allel.read_vcf(vcf_filename, samples=pop3)
    #sorted(callset.keys())
    gt_pop1 = allel.GenotypeArray(callset_pop1['calldata/GT'])
    gt_pop2 = allel.GenotypeArray(callset_pop2['calldata/GT'])
    gt_pop3 = allel.GenotypeArray(callset_pop3['calldata/GT'])
    #
    ac1 = gt_pop1.count_alleles()
    ac2 = gt_pop2.count_alleles()
    ac3 = gt_pop3.count_alleles()
    #
    pbs_test = allel.pbs(ac1, ac2, ac3, window_size=w_size, window_start=w_start, window_stop=w_stop, window_step=w_step, normed=True)
    #
    windows = list(allel.index_windows(np.arange(len(callset_pop1['variants/POS'])), size=w_size, start=w_start, stop=w_stop, step=w_step))
    window_chr = callset_pop1['variants/POS'][np.array(windows)-[0,1]]
    #
    chr = callset_pop1['variants/CHROM'][np.array(windows)[:,0]]
    #
    tag_chr = '_'.join(pathlib.Path(vcf_filename).name.split('_')[:2])
    tag_pop1 = re.findall(r"pop1_(.*?)\.",pop1_filename)[0]
    df=pd.DataFrame({"#chr":chr,"start":window_chr[:,0],"end":window_chr[:,1],"avg_pos":(window_chr[:,0]+window_chr[:,1])/2,"score":np.clip(pbs_test,0,None)})
    df["population"]=re.findall("pop1_([^_.]*).txt", pop1_filename)[0]
    df["test"]="pbs" # CHANGE THE THEST HERE
    df.to_csv(f"{tag_chr}_{tag_pop1}_pbs.tsv", sep='\t', index=False)
    #


file_pairs = list(itertools.product(vcf_files, pop1_files))

task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
# run_vcf_file(vcf_files[task_id])
run_vcf_file(*file_pairs[task_id])
