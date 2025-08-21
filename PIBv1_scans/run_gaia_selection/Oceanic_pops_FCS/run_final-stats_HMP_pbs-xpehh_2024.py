#!/usr/bin/env python3


import pandas as pd
import numpy as np
import scipy
import os
import pathlib
from glob import glob
from typing import Optional, Tuple

# CRATE THE WINDOWNS, WITH THE QUANTILES 


pbs = sorted(glob("../Oceanic_pops_PBS/*final.tsv"))
xpehh = sorted(glob("../Oceanic_pops_XPEHH/*final.tsv"))

pop_files = list(zip(pbs, xpehh))

class CustomIndexer(pd.api.indexers.BaseIndexer):
	def __init__(self, custom_start, custom_end, *args, **kwargs):
		super(CustomIndexer, self).__init__(*args, **kwargs)
		self._custom_start = custom_start
		self._custom_end = custom_end
	def get_window_bounds(#self, num_values, min_periods, center, closed, step
		self,
		num_values: int = 0,
		min_periods: Optional[int] = None,
		center: Optional[bool] = None,
		closed: Optional[str] = None,
		step: Optional[int] = None,
	) -> Tuple[np.ndarray, np.ndarray]:
		return self._custom_start, self._custom_end

#pbs_file, xpehh_file, = pop_files[0]

# for pbs_file, xpehh_file, in pop_files: 
def foo(pbs_file, xpehh_file,):
	pbs_df = pd.read_csv(pbs_file, sep='\t')
	xpehh_df = pd.read_csv(xpehh_file, sep='\t', usecols=['chr','avg_pos','end','quantile_rank'])
	## xpehh
	where_start = xpehh_df.reset_index().set_index(['chr','avg_pos']).loc[pbs_df[['chr','start']].to_numpy().tolist()]['index'].to_numpy()
	where_end = xpehh_df.reset_index().set_index(['chr','end']).loc[pbs_df[['chr','end']].to_numpy().tolist()]['index'].to_numpy()
	## max_scores = np.array([xpehh_df.loc[s_i:e_i+1].score.max() for s_i, e_i in np.vstack([where_start, where_end]).T])  ## too computationally heavy
	max_quantile_xpehh = xpehh_df.quantile_rank.rolling(CustomIndexer(xpehh_df.index.to_numpy(), np.repeat(where_end+1, np.hstack([np.diff(where_start), xpehh_df.shape[0]-where_start[-1]])))).max().loc[where_start]
	pbs_df['max_quantile_rank_xpehh'] = max_quantile_xpehh.to_numpy()
	#Harmonic mean p-val
	pbs_df['HMP'] = scipy.stats.hmean([1 - pbs_df['quantile_rank'], 1 - pbs_df['max_quantile_rank_xpehh']] , weights=None, nan_policy='omit')
	#pbs_df2 = pbs_df[pbs_df['HMP']<0.01]
	#
	new_file = pathlib.Path(pbs_file).name.replace("_final.tsv","_final_quantiles_window_stats_HMP_pbs-xpehh.tsv")
	#new_file2 = pathlib.Path(pbs_file).name.replace("_final.tsv","_final_quantiles_window_stats_HMP_pbs-xpehh_0.01.tsv")
	#
	pbs_df.to_csv(f"{new_file}", sep='\t', index=False)
	#pbs_df2.to_csv(f"{new_file2}", sep='\t', index=False)


task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
foo(*pop_files[task_id])
