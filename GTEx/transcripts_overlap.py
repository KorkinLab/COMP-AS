import numpy as np

def read_transcripts_list(fname):
	transcripts = []
	with open(fname) as f:
		for line in f:
			transcripts.append(line.split('.')[0])
	return transcripts[1:]

transcripts_gtex = read_transcripts_list('gtex-transcripts.txt')
transcripts_comp_as = read_transcripts_list('comp-as-transcripts.txt')

overlap = np.intersect1d(transcripts_gtex, transcripts_comp_as)
print(len(overlap))
print(len(overlap)/len(transcripts_comp_as))

with open('transcripts_overlap.txt', 'w') as f:
	for transcript in overlap:
		f.write('{}\n'.format(transcript))

def write_lost_transcripts(name, transcripts_list, overlap):
	with open('transcripts_lost_{}.txt'.format(name), 'w') as f:
		set_overlap = set(overlap)
		for transcript in transcripts_list:
			if transcript not in set_overlap:
				f.write('{}\n'.format(transcript))

write_lost_transcripts('comp_as', transcripts_comp_as, overlap)
write_lost_transcripts('gtex', transcripts_gtex, overlap)
