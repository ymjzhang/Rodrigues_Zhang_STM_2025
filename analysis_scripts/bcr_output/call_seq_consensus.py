# created by Duncan Morgan 
# modified by Jason Zhang
# last modified on Dec 21 2024
import pandas as pd
import numpy as np
import sys
#from Bio.Alphabet import generic_dna, Gapped
from Levenshtein import distance
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import Bio.Align
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio.Align.AlignInfo import SummaryInfo

[program, inputfile, outfile] = sys.argv
seq_input = pd.read_csv(inputfile, sep = '\t', index_col = 0)

seq_input['CONSENSUS_SEQUENCE'] = 'NA'
seq_input['CONSENSUS_SEQCOUNT'] = 'NA'
print(seq_input.head())
i = 0
for curr_bc in set(seq_input.LANE_ID):
    i += 1
    curr_data = seq_input.loc[seq_input.LANE_ID == curr_bc]
    curr_data = curr_data.sort_values(by = 'R2CONSCOUNT', ascending = False)
    consensus_seq = curr_data.SEQUENCE_IMGT.iloc[0]
    seqs_subset = curr_data
    if (curr_data.shape[0] > 1):
        
        
        X = pdist(np.array(curr_data.SEQUENCE_IMGT.tolist(), dtype=object).reshape(-1, 1), lambda x, y: distance(str(x[0]), str(y[0])))
	#X = pdist(np.array([x for x in curr_data.SEQUENCE_IMGT]).reshape(-1,1), lambda x, y: distance(x[0], y[0]))
        Z = hierarchy.linkage(X, 'single')
        clusters = hierarchy.cut_tree(Z, height = 5).flatten()

        seqs_subset = curr_data[clusters == Counter(clusters).most_common()[0][0]]
        seqs_subset.loc[:, 'LENGTH'] = [len(x) for x in seqs_subset.SEQUENCE_IMGT]
        commonLen = Counter(seqs_subset.LENGTH).most_common()[0][0]
        
        incomplete = True
        seqs_subset = seqs_subset[seqs_subset.LENGTH == commonLen].sort_values(by = 'R2CONSCOUNT', ascending = False)
        
        while(incomplete):
            seqs = list([SeqRecord(Seq(seqs_subset.SEQUENCE_IMGT.iloc[x])) for x in range(0,seqs_subset.shape[0])])
            seqs = Bio.Align.MultipleSeqAlignment(seqs)
            summary = SummaryInfo(seqs)
            consensus_seq = str(summary.gap_consensus(ambiguous = 'N', threshold = .5))
            
            if consensus_seq.count('N') > 0 & (seqs_subset.shape[0] > 1):
                seqs_subset = seqs_subset[:-1]
            else:
                incomplete = False

        # add handling for ns
        #if consensus_seq.count('N') > 0:
        #    break
    seq_input.loc[curr_data.index, 'CONSENSUS_SEQUENCE'] = consensus_seq
    seq_input.loc[curr_data.index, 'CONSENSUS_SEQCOUNT'] = seqs_subset.shape[0]

    if i % 100 == 0:
        print(i)
		
		
seq_input['ERRORDIST'] = [distance(seq_input.CONSENSUS_SEQUENCE.iloc[x], seq_input.SEQUENCE_IMGT.iloc[x]) for x in range(0, seq_input.shape[0])]
seq_input['CONSENSUS_AMBIG'] = [str(seq_input.CONSENSUS_SEQUENCE.iloc[x]).count('N') for x in range(0, seq_input.shape[0])]
seq_input.to_csv(outfile, sep = '\t')
