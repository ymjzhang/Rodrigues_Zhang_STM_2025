# created by Duncan Morgan 
# Last updated by Jason Zhang
# Last updated on 10/27/2022 

import pandas as pd
import numpy as np
import os

metafile = 'metadata.csv' #### Change line 
seq_project = '221004Lov' ##### Change Line 

os.system('mkdir BCLists')

metadata = pd.read_csv(metafile)
sample_id = pd.read_csv('hto_samples.csv', sep = ',')
samples = set(metadata.orig)
print(sample_id)
print(metadata.head())

for curr_sample in samples:
	print(curr_sample)
	curr_BC = metadata.bc.loc[metadata.orig == curr_sample].values
	curr_id = sample_id.DNAid.loc[sample_id.Sample == curr_sample].values[0]
	curr_bc_file = 'BCLists/' + curr_sample + '.txt'
	
	x = filter(lambda x:'N' not in x, curr_BC)
	np.savetxt(curr_bc_file, list(x), fmt = '%s')

	fastq1 = seq_project + '_' + curr_id + '_1_sequence.fastq.gz'
	fastq2 = seq_project + '_' + curr_id + '_2_sequence.fastq.gz'

	command_text = 'CITE-seq-Count -R1 ' + fastq1 + ' -R2 ' + fastq2 + ' -cbf 1 -cbl 12 -umif 13 -umil 20 -t hash_bcs_mouse.csv -cells 8000 -wl ' + curr_bc_file + ' -o Output/' + curr_sample + ' --dense' 
	print(command_text)
	os.system(command_text)