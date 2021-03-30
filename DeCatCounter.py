#!/usr/bin/env python

"""
This is DeCatCounter, a DeconCatenator and Counter for CCS sequencing data

Author: Celia Blanco 
Email: celiablanco@ucla.edu

V1.0 2021-03-01 : Created	

Dependencies:
	• python
	• biopython
	• python-Levenshtein
	• tabulate
	• pandas

How to run it:
$ python DeCatCounter.py in_file barcodes.txt tol_dem_f tol_dem_r constant.txt tol_decon_f tol_decon_r translation(y/n) len1 len2

e.g.:
$ python DeCatCounter.py CCS_99.fasta barcodes.txt 2 2 constant.txt 4 2 y 707 825 

"""

import sys
import os
import math
import time
import Levenshtein
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO
from Bio.Seq import Seq
from tabulate import tabulate


print("""
 ____        ____      _    ____                  _            
|  _ \  ___ / ___|__ _| |_ / ___|___  _   _ _ __ | |_ ___ _ __ 
| | | |/ _ \ |   / _` | __| |   / _ \| | | | '_ \| __/ _ \ '__|
| |_| |  __/ |__| (_| | |_| |__| (_) | |_| | | | | ||  __/ |   
|____/ \___|\____\__,_|\__|\____\___/ \__,_|_| |_|\__\___|_|                                                           

""")

# genetic code - dictionary used for translation
gencode = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
		
def splitAt(pos,split,lev):
	"""Split a sequence based on constant region positions"""
	poss = [1]
	for element in split:
		poss.append(element + 1)
	poss.append(len(pos) + 1)
	out_pos = []
	out_lev = []
	out_lev_min = []
	out_winner = []
	for i in range (0, len(poss)-1):
		out_pos.append(pos[poss[i]-1:poss[i+1]-1])
		out_lev.append(lev[poss[i]-1:poss[i+1]-1])
	for i, element in enumerate(out_lev):
		min_value = min(element)
		out_lev_min.append([j for j, x in enumerate(element) if x==min_value])
	for i, element in enumerate(out_pos):
		out_winner.append([])
		for element2 in out_lev_min[i]:
			out_winner[i].append(out_pos[i][element2])
	return(out_winner)

def translate_codon(codon):
	"""Translate a single codon"""
	return gencode.get(codon.upper(), 'x')
 
def split_into_codons(dna, frame):
	"""Split a sequence into codons"""
	codons = []
	for i in range(frame - 1, len(dna)-2, 3):
		codon = dna[i:i+3]
		codons.append(codon)
	return codons
 
def translate_dna_single(dna, frame=1):
	"""Translate a dna sequence in a single frame"""
	codons = split_into_codons(dna, frame)
	amino_acids = ''
	for codon in codons:
		if translate_codon(codon) == "*":
			amino_acids = amino_acids + translate_codon(codon)
			break
		else:
			amino_acids = amino_acids + translate_codon(codon)
	return amino_acids


def main():
	############################################################
	############### PARAMETERS DEMULTIPLEX ###############
	############################################################

	# define demultiplexing parameters
	file_in_to_dem = sys.argv[1]
	bc_file = sys.argv[2]
	tol_dem_f = int(sys.argv[3])
	tol_dem_r = int(sys.argv[4])

	# generate pair of barcodes assigned to each sample
	bcs = open(bc_file, 'r')
	samples = []
	bc_f = []
	bc_r = []
	bc_pair_fdir = []
	bc_pair_rdir = []
	for line in bcs:
		linesp = line.split()
		samples.append(linesp[0])
		bc_f.append(linesp[1])
		bc_r.append(linesp[2])
		bc_pair_fdir.append((linesp[1],linesp[2]))
		bc_pair_rdir.append((str(Seq(linesp[2]).reverse_complement()),str(Seq(linesp[1]).reverse_complement())))

	############################################################
	############### PARAMETERS DECONCATENATE ###############
	############################################################

	# function to translate a single codon
	ad_file = sys.argv[5]
	tol_dec_f = int(sys.argv[6])
	tol_dec_r = int(sys.argv[7])

	# assign constant regions for fwd and rev seq reads
	ads = open(ad_file, 'r')
	adapters=[]
	for line in ads:
		linesp = line.split()
		adapters.append(linesp[0])

	ad_f = adapters[0]
	ad_r = adapters[1]
	ad_f_RC = str(Seq(adapters[0]).reverse_complement())
	ad_r_RC = str(Seq(adapters[1]).reverse_complement())

	############################################################
	############### PARAMETERS FILTER ###############
	############################################################

	# set lenght thresholds (e.g. +/- 5% full length final sequence)
	low_len = int(sys.argv[9])
	hi_len = int(sys.argv[10])

	# read input file
	print("")
	recorcito_dem = []
	total_in = 0
	with open(file_in_to_dem, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		if any(fasta) == True:
			print("Reading input FASTA file ...")
			fastx_sequences = SeqIO.parse(open(file_in_to_dem),"fasta")
		else:
			handle.seek(0)
			fastq = SeqIO.parse(handle, "fastq")
			if any(fastq) == True:
				print("Reading input FASTQ file ...")
				fastx_sequences = SeqIO.parse(open(file_in_to_dem),"fastq")
			else:
				print("""Input file is neither FASTA nor FASTQ. No further action will be performed. Please provide input file in either FASTA or FASTQ format
				""")
				sys.exit(0)
	for fasta in fastx_sequences:
			name, sequence = fasta.id, str(fasta.seq)
			recorcito_dem.append({'name': name, 'seq': sequence})
			total_in += 1

	#generate list to stored demultiplexed sequences (will be the input for generating count files)
	dem_seq = []
	for sample in samples:
		dem_seq.append(0)
	
	#generate list to stored deconcatenated sequences (will be the input for generating count files)
	dec_seq = []
	for sample in samples:
		dec_seq.append(0)
		
	#generate list to stored trimmed sequences (will be the input for generating count files)
	trim_seq = []
	for sample in samples:
		trim_seq.append([])
	
	### demultiplex fwd reads
	print("Processing forward reads ...")
	print("")
	for element in recorcito_dem:
		for s,pair in enumerate(bc_pair_fdir):
			if Levenshtein.distance(element['seq'][0:len(pair[0])],pair[0]) <= tol_dem_f:
				if Levenshtein.distance(element['seq'][-len(pair[1]):],pair[1]) <= tol_dem_r:
					dem_seq[s] += 1
	# deconcatenate fwd reads
					seq_temp = element['seq'][len(pair[0]):len(element['seq'])-len(pair[1])]
					vec_seq_lev = []
					vec_seq_dis = []
					vec_pos = []
					vec_lev = []
					vec_dif = []
					vec_spl = []
					vec_loc = []
					vec_sta = []
					vec_end = []
					vec_adaseq = []
					vec_mat = []
					vec_ins = []
					for i in range(0, len(seq_temp)+1):
						lev = Levenshtein.distance(seq_temp[i:i+len(ad_f)],ad_f)
						vec_seq_lev.append(lev)
						if lev <= tol_dec_f:
							vec_seq_dis.append(int(1))
							vec_pos.append(i+1)
							vec_lev.append(lev)
						else:
							vec_seq_dis.append(0)
					if 1 in vec_seq_dis:
						for j in range(0, len(vec_pos)-1):
							vec_dif.append(vec_pos[j+1]-vec_pos[j])
						for k in range(0, len(vec_dif)):
							if vec_dif[k] > 1:
								vec_spl.append(k+1)
						spl = splitAt(vec_pos,vec_spl, vec_lev)
						for m in range(0, len(spl)):
							vec_loc.append(spl[m][int(math.floor(len(spl[m])/2))])
							vec_sta.append(spl[m][int(math.floor(len(spl[m])/2))])
						for startpos in vec_sta:
							vec_end.append(startpos + len(ad_f) - 1)
							vec_adaseq.append(seq_temp[startpos-1:startpos + len(ad_f) - 1])
						if vec_sta[0] != 1:
							vec_mat.append([1,vec_sta[0]-1])
						for n in range(0, len(vec_sta)-1):
							vec_mat.append([vec_end[n]+1, vec_sta[n+1]-1])
						vec_mat.append([vec_end[-1]+1,len(seq_temp)])
						num = len(vec_mat)
						for p in vec_mat:
							vec_ins.append(seq_temp[p[0]-1:p[1]])
							dec_seq[s]+=1
	# trim fwd reads													
							seq_temp2 = str(seq_temp[p[0]-1:p[1]])
							vec_seq_lev2 = []
							vec_seq_dis2 = []
							vec_pos2 = []
							vec_lev2 = []
							vec_dif2 = []
							vec_spl2 = []
							vec_loc2 = []
							vec_sta2 = []
							vec_end2 = []
							vec_adaseq2 = []
							vec_mat2 = []
							vec_ins2 = []
							for q in range(0, len(seq_temp2)+1):
								lev2 = Levenshtein.distance(seq_temp2[q:q+len(ad_r)],ad_r)
								vec_seq_lev2.append(lev2)
								if lev2 <= tol_dec_r:
									vec_seq_dis2.append(int(1))
									vec_pos2.append(q+1)
									vec_lev2.append(lev2)
								else:
									vec_seq_dis2.append(0)
							if 1 in vec_seq_dis2:
								for t in range(0, len(vec_pos2)-1):
									vec_dif2.append(vec_pos2[t+1]-vec_pos2[t])
								for r in range(0, len(vec_dif2)):
									if vec_dif2[r] > 1:
										vec_spl2.append(r+1)
								spl2 = splitAt(vec_pos2,vec_spl2, vec_lev2)
								for w in range(0, len(spl2)):
									vec_loc2.append(spl2[w][int(math.floor(len(spl2[w])/2))])
									vec_sta2.append(spl2[w][int(math.floor(len(spl2[w])/2))])
								for startpos in vec_sta2:
									vec_end2.append(startpos + len(ad_r) - 1)
									vec_adaseq2.append(seq_temp2[startpos-1:startpos + len(ad_r) - 1])
								if vec_sta2[0] != 1:
									vec_mat2.append([1,vec_sta2[0]-1])	
								for g in range(0, len(vec_sta2)-1):
									vec_mat2.append([vec_end2[g]+1, vec_sta2[g+1]-1])
								vec_mat2.append([vec_end2[-1]+1,len(seq_temp2)])
								num2 = len(vec_mat2)
								for y in vec_mat2:
									insertion = seq_temp2[y[0]-1:y[1]]
									vec_ins2.append(insertion)
								if low_len <= len(vec_ins2[0]) <= hi_len:
									trim_seq[s].append(vec_ins2[0])
	
	# demultiplex rev reads
	print("Processing reverse reads ...")
	print("")
	for element in recorcito_dem:
		for s,pair in enumerate(bc_pair_rdir):
			if Levenshtein.distance(element['seq'][0:len(pair[0])],pair[0]) <= tol_dem_r:
				if Levenshtein.distance(element['seq'][-len(pair[1]):],pair[1]) <= tol_dem_f:
					dem_seq[s] += 1
	# deconcatenate rev reads
					seq_temp = element['seq'][len(pair[0]):len(element['seq'])-len(pair[1])]
					vec_seq_lev = []
					vec_seq_dis = []
					vec_pos = []
					vec_lev = []
					vec_dif = []
					vec_spl = []
					vec_loc = []
					vec_sta = []
					vec_end = []
					vec_adaseq = []
					vec_mat = []
					vec_ins = []
					for i in range(0, len(seq_temp)+1):
						lev = Levenshtein.distance(seq_temp[i:i+len(ad_r_RC)],ad_r_RC)
						vec_seq_lev.append(lev)
						if lev <= tol_dec_r:
							vec_seq_dis.append(int(1))
							vec_pos.append(i+1)
							vec_lev.append(lev)
						else:
							vec_seq_dis.append(0)
					if 1 in vec_seq_dis:
						for j in range(0, len(vec_pos)-1):
							vec_dif.append(vec_pos[j+1]-vec_pos[j])
						for k in range(0, len(vec_dif)):
							if vec_dif[k] > 1:
								vec_spl.append(k+1)
						spl = splitAt(vec_pos,vec_spl, vec_lev)
						for m in range(0, len(spl)):
							vec_loc.append(spl[m][int(math.floor(len(spl[m])/2))])
							vec_sta.append(spl[m][int(math.floor(len(spl[m])/2))])
						for startpos in vec_sta:
							vec_end.append(startpos + len(ad_r_RC) - 1)
							vec_adaseq.append(seq_temp[startpos-1:startpos + len(ad_r_RC) - 1])
						if vec_sta[0] != 1:
							vec_mat.append([1,vec_sta[0]-1])
						for n in range(0, len(vec_sta)-1):
							vec_mat.append([vec_end[n]+1, vec_sta[n+1]-1])
						vec_mat.append([vec_end[-1]+1,len(seq_temp)])
						num = len(vec_mat)
						for p in vec_mat:
							vec_ins.append(seq_temp[p[0]-1:p[1]])
							dec_seq[s]+=1
	# trim rev reads
							seq_f = str(seq_temp[p[0]-1:p[1]])
							seq_r = str(Seq(seq_f).reverse_complement())
							seq_temp2 = str(seq_temp[p[0]-1:p[1]])
							vec_seq_lev2 = []
							vec_seq_dis2 = []
							vec_pos2 = []
							vec_lev2 = []
							vec_dif2 = []
							vec_spl2 = []
							vec_loc2 = []
							vec_sta2 = []
							vec_end2 = []
							vec_adaseq2 = []
							vec_mat2 = []
							vec_ins2 = []
							for q in range(0, len(seq_temp2)+1):
								lev2 = Levenshtein.distance(seq_temp2[q:q+len(ad_f_RC)],ad_f_RC)
								vec_seq_lev2.append(lev2)
								if lev2 <= tol_dec_f:
									vec_seq_dis2.append(int(1))
									vec_pos2.append(q+1)
									vec_lev2.append(lev2)
								else:
									vec_seq_dis2.append(0)
							if 1 in vec_seq_dis2:
								for t in range(0, len(vec_pos2)-1):
									vec_dif2.append(vec_pos2[t+1]-vec_pos2[t])
								for r in range(0, len(vec_dif2)):
									if vec_dif2[r] > 1:
										vec_spl2.append(r+1)
								spl2 = splitAt(vec_pos2,vec_spl2, vec_lev2)
								for w in range(0, len(spl2)):
									vec_loc2.append(spl2[w][int(math.floor(len(spl2[w])/2))])
									vec_sta2.append(spl2[w][int(math.floor(len(spl2[w])/2))])
								for startpos in vec_sta2:
									vec_end2.append(startpos + len(ad_f_RC) - 1)
									vec_adaseq2.append(seq_temp2[startpos-1:startpos + len(ad_f_RC) - 1])
								if vec_sta2[0] != 1:
									vec_mat2.append([1,vec_sta2[0]-1])
								for g in range(0, len(vec_sta2)-1):
									vec_mat2.append([vec_end2[g]+1, vec_sta2[g+1]-1])
								vec_mat2.append([vec_end2[-1]+1,len(seq_temp2)])
								num2 = len(vec_mat2)
								for y in vec_mat2:
									insertion = seq_temp2[y[0]-1:y[1]]
									vec_ins2.append(insertion)
								if low_len <= len(vec_ins2[0]) <= hi_len:
									seq_f = vec_ins2[0]
									seq_r = str(Seq(seq_f).reverse_complement())
									trim_seq[s].append(seq_r)
	
	# create directory for count files
	print("Creating output directory ...")
	print("")
	timestr = time.strftime("%Y%m%d-%H%M%S")
	if not os.path.exists("output" + str(timestr) + "/counts"):
		os.makedirs("output" + str(timestr) + "/counts")
	
	# create count files	
	f_name_out_counts = []
	for sample in samples:
		f_name_out_counts.append("output" + str(timestr) + "/counts"
		+"/sample" + str(sample) + ".txt")

	# open count files to write
	out_counts = []
	for file in f_name_out_counts:
		out_counts.append(open(file, 'w'))

	# read translation flag
	trans_flag = str(sys.argv[8])

	# create directory for aa count files
	if trans_flag == 'y':
		if not os.path.exists("output" + str(timestr) + "/counts_aa"):
			os.makedirs("output" + str(timestr) + "/counts_aa")

	# create aa count files	
		f_name_out_counts_aa = []
		for sample in samples:
			f_name_out_counts_aa.append("output" + str(timestr) + "/counts_aa"
			+"/sample" + str(sample) + "_aa" + ".txt")

	# open aa count files
		out_counts_aa = []
		for file in f_name_out_counts_aa:
			out_counts_aa.append(open(file, 'w'))

	# generate content for count files
	for i, sample in enumerate(samples):
		print("Generating count file for sample " + str(sample) + " ...")
		print("")
		df = DataFrame (trim_seq[i],columns=['seq'])
		df = df.groupby('seq').seq.count().reset_index(name='counts')
		df = df.sort_values(by=['counts'], ascending=False)
		df['aa'] = df['seq'].apply(translate_dna_single)
		df_nt = df[['seq', 'counts']].copy()
		nmol = df_nt['counts'].sum()
		nseq = len(df_nt.index)
		out_counts[i].write("number of unique sequences = "+ str(nseq)+ '\n')
		out_counts[i].write("total number of molecules = "+ str(nmol)+ '\n')
		out_counts[i].write('\n')
		out_counts[i].write(df_nt.to_csv(header=False, index=False, sep="\t"))

	# translate to amino acids	
		if trans_flag == 'y':
			print("Generating aa count file for sample " + str(sample) + " ...")
			print("")
			df_aa_dup = df[['counts', 'aa']].copy()
			cols_aa = ['aa', 'counts']
			df_aa_dup = df_aa_dup[cols_aa]
			df_aa = df_aa_dup.groupby('aa').counts.sum().reset_index(name='counts')
			df_aa = df_aa.sort_values(by=['counts'], ascending=False)
			nmol = df_aa['counts'].sum()
			nseq = len(df_aa.index)
			out_counts_aa[i].write("number of unique sequences = "+ str(nseq)+ '\n')
			out_counts_aa[i].write("total number of molecules = "+ str(nmol)+ '\n')
			out_counts_aa[i].write('\n')
			out_counts_aa[i].write(df_aa.to_csv(header=False, index=False, sep="\t"))

	# close count files
	for i, file in enumerate(f_name_out_counts):
		out_counts[i].close()

	# close amino acid count files
	if trans_flag == 'y':
		for i, file in enumerate(f_name_out_counts_aa):
			out_counts_aa[i].close()

	# create log ouput file
	f_out=open("output" + str(timestr) + "/log.txt", 'w') 

	f_out.write("Input file: " + str(file_in_to_dem)+ '\n')	
	f_out.write("Barcodes file: " + str(bc_file)+ '\n')	
	f_out.write("Tolerance for forward barcode: "+ str(tol_dem_f)+ '\n')	
	f_out.write("Tolerance for reverse barcode: "+ str(tol_dem_r)+ '\n')	
	f_out.write("Constant regions file: "+ str(ad_file)	+ '\n')
	f_out.write("Tolerance for forward constant region: "+ str(tol_dec_f)+ '\n')	
	f_out.write("Tolerance for reverse constant region: "+ str(tol_dec_r)+ '\n')
	f_out.write("Lower length: "+ str(low_len)+ '\n')
	f_out.write("Higher length: "+ str(hi_len)+ '\n')
	f_out.write(""+ '\n')
	f_out.write("Number of sequences in input file: "+ str(total_in)+ '\n')	
	f_out.write(""+ '\n')	

	header = ["Sample", "Demultiplexed", "Deconcatenated", "Filtered"]
	table = []
	for i, sample in enumerate(samples):
		table.append([])
		table[i].append(samples[i])
		table[i].append(dem_seq[i])
		table[i].append(dec_seq[i])
		table[i].append(len(trim_seq[i]))

	print(tabulate(table, headers = header))
	f_out.write(tabulate(table, headers = header))
	f_out.close()
	
if __name__=='__main__':
	main()
