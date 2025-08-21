#!/usr/bin/env python3

import numpy as np
import pandas as pd
import re
import sys
import os
from glob import glob

files_genecode = sorted(glob("*closest.tsv"))

# for file in files_genecode:
def name_of_function(file):
        print(file)
        sys.stdout.flush()
        df=pd.read_csv(file , sep='\t', header=None)
        df2=df.rename(columns={0: "chr", 1: "start", 2: "end", 3: "avg_pos",4: "score", 5:"population", 6:"test", 10:"ensembl", 11:"gene", 12:"distance"})
        ## remove chr from chr column
        #df2["chr"]=df2["chr"].str.replace("chr", "").astype('int') # this line works if chr1 chr2 ...
        ####
        df3=df2[["chr","start","end","avg_pos","score","population","test","ensembl","gene","distance"]].reset_index()
        ## dealing with topgene below
        ###OLD LINE in pandas v1.x.x: temp_df = df3.reset_index()[np.isfinite(df3.score)].sort_values('score', ascending=False)
        temp_df = df3[np.isfinite(df3.score)].sort_values('score', ascending=False)
        df3['topgene'] = ''
        ###OLD LINE in pandas v1.x.x: top_genes = temp_df.groupby('gene').index.nth(0).reset_index()
        top_genes = temp_df.groupby('gene').gene.nth(0).reset_index()
        ###OLD LINE in pandas v2.x.x: df3.topgene.loc[top_genes['index'].to_numpy()] = top_genes['gene'].to_list()
        df3.loc[top_genes['index'].to_numpy(), 'topgene'] = top_genes['gene'].to_list()
        df3['is_topgene'] = df3.topgene.map(bool)
        ####
        df3_groups=df3.groupby(['chr','end'])
        ###OLD LINE in pandas v1.x.x: df4=df3_groups.nth(0).drop('is_topgene', axis=1)
        df4=df3_groups.nth(0).drop('is_topgene', axis=1).set_index(['chr','end'])
        df4.gene = df3_groups.gene.apply(lambda L: ",".join(L)) # ***
        foo = df3_groups.is_topgene.sum().map(bool)
        ###OLD LINE in pandas v1.x.x: df4.topgene[foo] = df3[df3.is_topgene].groupby(['chr','end']).topgene.apply(lambda L: ",".join(L[L!=""]))
        df4.loc[foo, 'topgene'] = df3[df3.is_topgene].groupby(['chr','end']).topgene.apply(lambda L: ",".join(L[L!=""]))
        df4=df4.reset_index().set_index("index").sort_index()
        df5=df4.drop("start",axis=1)
        df5.insert(1,"start",df4.start)
        df5["quantile_rank"] = df5["score"].rank(pct = True)
        new_file=file.replace('.tsv','')
        df5.to_csv(f"{new_file}_final.tsv", sep='\t', index=False)


task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))

name_of_function(files_genecode[task_id])
