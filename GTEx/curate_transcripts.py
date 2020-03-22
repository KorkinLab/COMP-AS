import pandas as pd
import numpy as np

dbs = ['tissues_SMTSD_tpm_mean.tsv', 'tissues_SMTS_tpm_mean.tsv']

transcripts = pd.read_csv('transcripts_overlap.txt', header=None)
transcripts = list(transcripts[transcripts.columns[0]])

for db in dbs:
	name = db.split('.')[0]
	data = pd.read_csv(db, sep='\t', index_col=0)
	#out_data = data.loc[transcripts,:]
	#out_data.to_csv('{}_curated.tsv'.format(name))
	#print(out_data)
	with open('duplicates.txt', 'w') as f:
		for x in data.index[data.index.duplicated()]:
			f.write('{}\n'.format(x))
	print(len(np.intersect1d(transcripts, data.index)))
