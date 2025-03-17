import pandas as pd
import numpy as np
import glob

files = glob.glob('Output/*/dense*')

for curr in files:
	print(curr)
	sample_name = curr.split('/')[1]
	
	curr_df = pd.read_csv(curr, sep = '\t', index_col = 0)
	curr_df.columns = sample_name + '_' + curr_df.columns
	if curr == files[0]:
		df = curr_df
	
	else:
		df = df.join(curr_df)
		
df = df.transpose()
df.to_csv('citeseq_results_dense.csv')