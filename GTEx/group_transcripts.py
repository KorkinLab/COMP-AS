import pandas as pd
import numpy as np

attr = pd.read_csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')
reader = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct', sep='\t', chunksize=10000)
chunks = []
for chunk in reader:
	chunks.append(chunk)
expressions = pd.concat(chunks)


sample_names = attr['SAMPID']
terms = ['SMTS', 'SMTSD']


for term in terms:
	samples_sum = 0
	groups = np.unique(attr[term])
	df = pd.DataFrame()
	for group in groups:
		group_samples = sample_names.loc[attr[term] == group]
		group_samples = [x.rstrip() for x in group_samples]
		#print(group)
		#print(group_samples)
		#intersect = group_samples
		intersect = np.intersect1d(group_samples, list(expressions.columns))
		#print(len(intersect)/len(group_samples))
		samples_sum += len(intersect)
		df[group] = np.mean(expressions[intersect], axis=1)
	print('{}: {} samples'.format(term, samples_sum))
	index=[x.split('.')[0] for x in expressions['transcript_id']]
	df.index = index #expressions['transcript_id']
	df.to_csv('tissues_{}_tpm_mean.tsv'.format(term), sep='\t')

