import os
import glob
files = glob.glob('*R2*.fastq.gz') + glob.glob('*R3*.fastq.gz')

for curr in files:
	command = 'sbatch run_fastqc.sh ' + curr
	print(command)
	os.system(command)
