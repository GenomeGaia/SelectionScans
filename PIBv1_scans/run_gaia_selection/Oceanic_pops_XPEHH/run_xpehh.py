#!/usr/bin/env python3

import os
import re
import pathlib
import itertools
import allel
import numpy as np
import pandas as pd
from glob import glob
from numpy import asarray
from numpy import savetxt

map_files = sorted(glob("../../key_files/gmap_files/chr*.b37.gmap"),key=lambda name: int(re.findall(r"chr([0-9]*?)[_\.]",name)[0]))

pop1_files = sorted(glob("pops_oceania/pop1_*.txt"))

file_name = "pops_oceania/pop2_EAS.txt"
with open(file_name,'r') as f:
    pop2 = f.read().split("\n")[:-1]

#vcf_filename = 'PIBv1_chr21_phased.MAF0.01_SNP_only.vcf.gz' # if True: (and remove the for loop)
vcf_files = sorted(glob('../../vcf_files/PIBv1_chr*_final_phased.SNP_only.vcf.gz'),key=lambda name: int(re.findall(r"chr([0-9]*?)[_\.]",name)[0]))

chr_file_pairs = zip(map_files,vcf_files)

# for vcf_filename in vcf_files:
def run_vcf_file(vcf_filename, map_filename, pop1_filename):
    with open(pop1_filename,'r') as f:
        pop1 = f.read().split("\n")[:-1]
    callset_pop1 = allel.read_vcf(vcf_filename, samples=pop1)
    callset_pop2 = allel.read_vcf(vcf_filename, samples=pop2)
    #Genotype Array
    gt_pop1 = allel.GenotypeArray(callset_pop1['calldata/GT'])
    gt_pop2 = allel.GenotypeArray(callset_pop2['calldata/GT'])
    # n_variants, n_haplotypes
    gt_h1 = gt_pop1.haploidify_samples()
    gt_h2 = gt_pop2.haploidify_samples()
    # haplotypes array 
    h1 = allel.HaplotypeArray(gt_h1, dtype='i1')
    h2 = allel.HaplotypeArray(gt_h2, dtype='i1')
    # 
    chr = callset_pop1['variants/CHROM']
    pos = callset_pop1['variants/POS']
    #
    # linear interpolation | genetic map
    map_data = pd.read_csv(map_filename, delimiter='\t')
    xp = map_data["pos"] # [0]
    fp = map_data["cM"] #[2]
    map_pos = np.interp(pos, xp, fp)
    #
    #XPEHH test
    #with np.errstate(divide='ignore'):
    xp = allel.xpehh(h1, h2, pos, map_pos=map_pos, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
    #
    xp[~np.isfinite(xp)] = np.nan    
    #
    score = allel.standardize(xp)
    #
    tag_chr = '_'.join(pathlib.Path(vcf_filename).name.split('_')[:2])
    tag_pop1 = re.findall(r"pop1_(.*?)\.",pop1_filename)[0]
    df=pd.DataFrame({"#chr":chr,"start":pos-1,"end":pos,"avg_pos":pos,"score":score})
    df["population"]=re.findall("pop1_([^_.]*).txt", pop1_filename)[0]
    df["test"]="xpehh" # CHANGE THE THEST HERE
    df.to_csv(f"{tag_chr}_{tag_pop1}_xpehh.tsv", sep='\t', index=False)
    #

all_files = list(itertools.product(list(zip(vcf_files, map_files))[:], pop1_files))
file_triplets = [(a,b,c) for (a,b),c in all_files]

task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
#run_vcf_file(vcf_files[task_id])
run_vcf_file(*file_triplets[task_id])
