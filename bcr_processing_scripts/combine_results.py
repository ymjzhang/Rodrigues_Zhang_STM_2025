import pandas as pd
import numpy as np
import glob

files = glob.glob('*/*final_call*')
files.sort()
for curr in files:
	print(curr)
	sample_name = curr.split('/')[0]
	curr_df = pd.read_csv(curr, sep = '\t')
	curr_df['SAMPLE'] = sample_name
	curr_df['ID'] = [curr_df.SAMPLE[x] +"_"+ curr_df.SEQUENCE_ID[x] for x in range(0, curr_df.shape[0])]
	curr_df.index = curr_df.ID
	curr_df = curr_df.loc[curr_df.R2CONSCOUNT>1,:]

	
	if curr == files[0]:
		df = curr_df
	else:
		df = pd.concat([df, curr_df])

df.to_csv('combined_results.tab', sep = '\t', index = False)
