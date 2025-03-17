import os
import glob

files = glob.glob('*/*.fastq')
for curr in files:
	command = 'sbatch zip.sh ' + curr
	print(command)
	os.system(command)
