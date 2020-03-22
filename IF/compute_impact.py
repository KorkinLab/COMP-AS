import numpy as np
import pandas as pd
from scipy.special import expit


def class ImpactFactor:
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
			'gene-max': self.max_gene
			'gene-rearrangement': self.rearrangement_gene
		}

	def sigmoid(x):
		return expit(x)

	# Implementation of "Use" metric from Fig.5
	def use(tissues_expressed_num, avg_transcript_rank, transcripts_num):
		return (tissues_expressed_num * avg_transcript_rank) / transcripts_num

	# Implementation of inverted "Extent rearrangement" metric from Fig.5
	def extent_rearrangement(modification_extent):
		bs_total = len(modification_extent)
		return bs_total / (bs_total - sigmoid(np.ones(bs_total) - np.array(modification_extent)))
	

	#########
	# Wrappers that execute IF function based on passed strategy
	# transcript_expr: expression level of the transcript of interest (TOI). Appropriate value should be used for tissue and transcript levels
	# gene_expr: expression level of all transcripts from gene of interest (GOI)
	# avg_transcript_rank: average rank of the transcript of interest (TOI). Appropriate value should be used for tissue and transcript/gene levels
	# transcripts_num: number of transcripts derived from the same gene as the TOI
	# tissues_expressed_num: number of tissues where TOI is expressed
	# bs_modification_percentage: a list of modification level for each binding site in TOI
	def transcript_in_tissiue_if(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage, strategy='average'):
		func = self.implementations.get('tissue-{}'.format(strategy), lambda: "Invalid strategy")
		return func(transcript_expr, avg_transcript_rank, bs_total, bs_modification_percentage)

	def transcript_if(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		func = self.implementations.get('transcript-{}'.format(strategy), lambda: "Invalid strategy")
		return func(transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage)

	def gene_if(self, gene_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		func = self.implementations.get('gene-{}'.format(strategy), lambda: "Invalid strategy")
	##########

	# Averaging strategy
	def avg_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return 

	def avg_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		pass

	def avg_gene(self, gene_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		pass

	# Max strategy
	def max_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		pass

	def max_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		pass

	def max_gene(self, gene_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		pass

	# Rearrangement strategy
	def rearrangement_transcript_in_tissue(self, transcript_expr, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return (avg_transcript_rank/transcripts_num) * extent_rearrangement(bs_modification_percentage)

	def rearrangement_transcript(self, transcript_expr, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return use(tissues_expressed_num, avg_transcript_rank, transcripts_num) * extent_rearrangement(bs_modification_percentage)

	def rearrangement_gene(self, gene_expr, avg_transcript_rank, tissues_expressed_num, avg_transcript_rank, transcripts_num, bs_modification_percentage):
		return use(tissues_expressed_num, avg_transcript_rank, transcripts_num) * extent_rearrangement(bs_modification_percentage)


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

data = load_data('../../../GTEx/tissues_SMTSD_tpm_mean_curated.tsv', 1000)
print(data)
