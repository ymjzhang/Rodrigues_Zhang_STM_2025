import os
import glob

files = glob.glob('*.fastq.gz')
for curr in files:
	command = 'sbatch unzip.sh ' + curr + ' --exclude=c[1-20]'
	os.system(command)

