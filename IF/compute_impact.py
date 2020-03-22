import numpy as np
import pandas as pd
from scipy.special import expit


def sigmoid(x):
	return expit(x)

class ImpactFactor:
	strategy = None
	implementations = None

	def __init__(self, strategy='average'):
		self.strategy = strategy
		self.implementations = {
			'tissue-average': self.avg_transcript_in_tissue,
			'tissue-max': self.max_transcript_in_tissue,
			'tissue-rearrangement': self.rearrangement_transcript_in_tissue,
			'transcript-average': self.avg_transcript,
			'transcript-max': self.max_transcript,
			'transcript-rearrangement': self.rearrangement_transcript,
			'gene-average': self.avg_gene,
			'gene-max': self.max_gene,
			'gene-rearrangement': self.rearrangement_gene
		}


	# Implementation of "Use" metric from Fig.5
	def use(self, tissues_expressed_num, avg_transcript_rank, transcripts_num):
		return (tissues_expressed_num * avg_transcript_rank) / transcripts_num

	# Implementation of inverted "Extent rearrangement" metric from Fig.5
	def extent_rearrangement(self, modification_extent):
		bs_total = len(modification_extent)
		return bs_total / (bs_total - sum(sigmoid(np.ones(bs_total) - np.array(modification_extent))))
	

	#########
	# Wrappers that execute IF function based on passed strategy
	# transcript_expr: expression level of the transcript of interest (TOI). Appropriate value should be used for tissue and transcript levels
	# gene_expr: expression level of all transcripts from gene of interest (GOI)
	# avg_transcript_rank: average rank of the transcript of interest (TOI). Appropriate value should be used for tissue and transcript/gene levels
	# transcripts_num: number of transcripts derived from the same gene as the TOI
	# tissues_expressed_num: number of tissues where TOI is expressed
	# bs_modification_percentage: a list of modification level for each binding site in TOI. For gene level modified binding sites for all transcripts should be listed
	def transcript_in_tissiue_if(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		func = self.implementations.get('tissue-{}'.format(self.strategy), lambda: "Invalid strategy")
		return func(transcript_expr, avg_transcript_rank, transcripts_num,  bs_modification_percentage)

	def transcript_if(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		func = self.implementations.get('transcript-{}'.format(self.strategy), lambda: "Invalid strategy")
		return func(transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage)

	def gene_if(self, gene_expr, tissues_expressed_num, transcripts_num, bs_modification_percentage):
		func = self.implementations.get('gene-{}'.format(self.strategy), lambda: "Invalid strategy")
		return func(gene_expr, tissues_expressed_num, transcripts_num, bs_modification_percentage)
	##########

	# Averaging strategy
	def avg_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return transcript_expr * np.mean(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	def avg_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return transcript_expr * np.mean(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	def avg_gene(self, gene_expr, tissues_expressed_num, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return gene_expr * np.mean(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	# Max strategy
	def max_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return transcript_expr * np.max(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	def max_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return transcript_expr * np.max(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	def max_gene(self, gene_expr, tissues_expressed_num, transcripts_num, bs_modification_percentage):
		bs_total = len(bs_modification_percentage)
		return gene_expr * np.max(sigmoid(np.ones(bs_total) - np.array(bs_modification_percentage)))

	# Rearrangement strategy
	def rearrangement_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return (avg_transcript_rank/transcripts_num) * self.extent_rearrangement(bs_modification_percentage)

	def rearrangement_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return self.use(tissues_expressed_num, avg_transcript_rank, transcripts_num) * self.extent_rearrangement(bs_modification_percentage)

	def rearrangement_gene(self, gene_expr, tissues_expressed_num, transcripts_num, bs_modification_percentage):
		return (tissues_expressed_num / transcripts_num) *self.extent_rearrangement(bs_modification_percentage)


def load_data(path, chunk=None, sep=','):
	if chunk is None:
		return pd.read_csv(path, sep=sep)
	else:
		reader = pd.read_csv(path, sep=sep, chunksize=chunk)
		chunks = []
		for chunk in reader:
			chunks.append(chunk)
			break
		expressions = pd.concat(chunks)
		return expressions


imp = ImpactFactor(strategy='rearrangement')
tissue = imp.transcript_in_tissiue_if(10, 3, 5, [0, 0.3, 0.5, 0.7])
transcript = imp.transcript_if(15, 7, 3, 5, [0, 0.3, 0.5, 0.7])
gene = imp.gene_if(30, 7, 5, [0, 0.3, 0.5, 0.7])

print('Tissue: {}\nTranscript: {}\n Gene: {}'.format(tissue, transcript, gene))

#data = load_data('../../../GTEx/tissues_SMTSD_tpm_mean_curated.tsv', 1000)
