#!/usr/bin/python3
import os
from Bio import SeqIO
import statistics
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

def parse_3UTRs(infasta_path, outfasta_path, id_delimiter="-"):
	records_3UTRs = []
	for record in SeqIO.parse(infasta_path, "fasta"):
		id_split = record.id.split(id_delimiter)
		transcript_id = id_split[0]
		polyA_count = int(id_split[1])
		stop_codon = record.seq[:3]
		if stop_codon == "TAA":
			seq_3UTR = record.seq[3:-polyA_count]
			record.name = transcript_id
			record.seq = seq_3UTR
			record.description = ""
			records_3UTRs.append(record)
	SeqIO.write(records_3UTRs, outfasta_path, "fasta")
	return records_3UTRs

def count_codon(seq, codon, frame=0, seq_len=0):
	codon_count = 0
	if seq_len == 0:
		seq_len = len(seq)
	for i in range(frame,seq_len,3):
		triplet = seq[i:i+3]
		if triplet == codon:
			codon_count+=1
	return codon_count

def count_codon_all_frames(records, seq_len=0, codon="TAA"):
	codon_count_frame0 = 0
	codon_count_frame1 = 0
	codon_count_frame2 = 0
	for record in records:
		seq = record.seq
		codon_count_frame0 += count_codon(seq, codon, seq_len=seq_len, frame=0)
		codon_count_frame1 += count_codon(seq, codon, seq_len=seq_len, frame=1)
		codon_count_frame2 += count_codon(seq, codon, seq_len=seq_len, frame=2)
	codon_tuple = (codon_count_frame0, codon_count_frame1, codon_count_frame2)
	return codon_tuple

def prepare_statistics_old(records):
	codons = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT']
	seq_lens = [3, 6, 9, 12, 15, 0]
	for seq_len in seq_lens:
		for codon in codons:
			codon_tuple = count_codon_all_frames(records, seq_len=seq_len, codon=codon)
			print ("codon", codon, "first", seq_len, "nt of 3'UTRs", codon_tuple)

def prepare_len_statistics(fasta_path):
	len_statistics = {}
	lengths = []
	for record in SeqIO.parse(fasta_path, "fasta"):
		length = len(record)
		lengths.append(length)
		if length in len_statistics.keys():
			len_statistics[length] += 1
		else:
			len_statistics[length] = 1
	max_len = max(lengths)
	print("median length", statistics.median(lengths))
	print("average length", statistics.mean(lengths))
	return max_len

def prepare_empty_codon_dict(max_len):
	codon_dict = {}
	codons = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT']
	for codon in codons:
		codon_dict[codon] = {}
		for codon_n in range(0, max_len//3 + 1):
			codon_dict[codon][codon_n] = {}
			for frame in [0, 1, 2]:
				codon_dict[codon][codon_n][frame] = 0
	return codon_dict

def add_codon_statistics_record(record, codon_dict):
	seq = record.seq
	seq_len = len(seq)
	length = seq_len
	for frame in [0,1,2]:
		codon_n = 0
		for i in range(frame,length,3):
			codon_n += 1
			triplet = seq[i:i+3]
			if len(triplet) == 3:
				codon_dict[triplet][codon_n][frame] += 1
	return codon_dict

def get_all_codon_statistics(fasta_path):
	max_len = prepare_len_statistics(fasta_path)
	codon_dict = prepare_empty_codon_dict(max_len)
	for record in SeqIO.parse(fasta_path, "fasta"):
		codon_dict = add_codon_statistics_record(record, codon_dict)
	return codon_dict

def make_graph(codon_dataframe, codon, outdir, x_label):
	outpath = outdir + codon + ".png"
	palette = sns.color_palette("hls", 3)
	myplot = sns.lineplot(data=codon_dataframe,palette=palette, dashes=None, linewidth=1)
	plt.xlabel(x_label)
	plt.ylabel("Codon count")
	# graph_title = f"Summary distribution of {codon} codon in 3'UTRs"
	# plt.title(graph_title)
	plt.legend(labels=["frame 1","frame 2","frame 3"])
	plt.savefig(outpath, dpi=300)
	plt.close()
	return 0


def make_graphs(codon_dict, outdir, upto_nts=None, upto_triplets=None):
	for codon in codon_dict:
		codon_data = codon_dict[codon]
		if upto_nts:
			codon_data_nt_postitions = {}
			for triplet_position in codon_data:
				nt_position = triplet_position * 3
				codon_data_nt_postitions[nt_position] =  codon_data[triplet_position]
			codon_dataframe = pd.DataFrame.from_dict(codon_data_nt_postitions, orient='index')
			codon_dataframe = codon_dataframe.truncate(after=upto_nts)
			x_label = ("Position, nt")
		elif upto_triplets:
			for triplet_position in codon_data:
				for frame in (2,1,0):
					codon_data[triplet_position][frame+1] = codon_data[triplet_position].pop(frame)
			codon_dataframe = pd.DataFrame.from_dict(codon_data, orient='index')
			codon_dataframe = codon_dataframe.truncate(after=upto_triplets)
			x_label = ("Position, triplets")
		else:
			print ("No x unit defined!")
			return 0
		print (codon)
		make_graph(codon_dataframe, codon, outdir, x_label)
	return 0

def convert_to_df(codon_dict, max_triplets=10):
	print("making data frame")
	df_dict = {}
	i = 0
	for codon in codon_dict:
		codon_data = codon_dict[codon]
		for position in codon_data:
			if position <= max_triplets:
				position_data = codon_data[position]
				for frame in position_data:
					current_dict = {}
					current_dict["codon"] = codon
					current_dict["frame"] = frame + 1
					current_dict["position, triplets"] = position
					current_dict["count"] = position_data[frame]
					df_dict[i] = current_dict
					i+=1
	codon_df = pd.DataFrame.from_dict(df_dict, orient='index')
	return codon_df

def make_facet_grid_graph(codon_df, outpath, max_x=30):
	print("making facet grid graph")
	sns.set_theme(style="ticks")
	palette = sns.color_palette("hls", 3)
	grid = sns.FacetGrid(codon_df, col_wrap=4, col="codon", hue="frame", palette=palette, sharex=False, legend_out=False)
	x_step = int(max_x / 5)
	grid.set(xticks=range(0, max_x+1, x_step))
	grid.map(sns.lineplot, "position, triplets","count")
	grid.add_legend()
	plt.savefig(outpath, dpi=300)
	return 0

fasta_path = "/Users/vl18625/work/blasto_local/3_UTR/extracted_3UTRs.fasta"
outdir = "/Users/vl18625/work/blasto_local/3_UTR/codon_usage_3UTRs_figures/"
outpath = "/Users/vl18625/work/blasto_local/3_UTR/codon_usage_facegrid.png"
print ("Preparing codon statistics")
codon_dict = get_all_codon_statistics(fasta_path)

max_triplets = 100
make_graphs(codon_dict, outdir, upto_triplets=max_triplets)
# codon_df = convert_to_df(codon_dict,max_triplets=max_triplets)
# make_facet_grid_graph(codon_df, outpath, max_x=max_triplets)

