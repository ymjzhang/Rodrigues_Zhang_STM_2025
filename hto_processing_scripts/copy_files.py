import pandas as pd
import glob
import os

hto_samples = pd.read_csv('hto_samples.csv')

files = glob.glob('../*/*.fastq.gz')
for curr in files:
	if any(dna_id in curr for dna_id in hto_samples['DNAid']):
		command = 'sbatch copy_files.sh ' + curr 
		print(command)
		os.system(command)
