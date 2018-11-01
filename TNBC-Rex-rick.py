#!/usr/bin/env python3

from datetime import datetime
from pathlib import Path

import numpy as np
import getpass
import multiprocessing
import pandas as pd
import collections
import subprocess as sp
import argparse
import math
import gzip
import glob
import sys
import os

def bcl2fastq(run_dir):
	sp.run([
		'nohup', 
		'bcl2fastq', 
		'--runfolder-dir', 
		run_dir, 
		'--barcode-mismatches=0', 
		'--no-lane-splitting'
		])

def cufflinks(bam_dir, out_dir):
	"""
	bam_dir  - path to directory containing bam files
	out_dir  - path to desired output directory (fpkm files)
	n_proc   - number of processors for multithreading (string)
	gtf_path - gene annotation file
	"""
	bam_list = list(Path(bam_dir).glob('*.bam'))
	for bam in bam_list:
		name = str(Path(str(bam).split('.')[0]).name)
		outfile = str(Path(out_dir)/'_'.join([name, 'cufflinks']))
		sp.run([
			'cufflinks', 
			'-p', 
			n_proc, 
			'--library-type', 
			'fr-firststrand',
			'-G', 
			gtf_path, 
			'-o', 
			outfile, 
			str(bam)
			])

def cufflinks_merge(fpkm_dir, out_dir, run_name):
	"""
	Produces a log table (csv) from a directory containing 
	fpkm tracking files.
	"""
	out_dir = Path(out_dir)
	fpkm_list = list(Path(fpkm_dir).glob('**/*genes.fpkm_tracking'))
	final_df = pd.DataFrame()
	for file in fpkm_list:
		## Change this
	    name = Path(file).parts[-2]
	    df = pd.read_table(file)
	    df = df[['gene_short_name', 'FPKM']].groupby('gene_short_name').mean()
	    df.rename(columns={'FPKM': name}, inplace=True)
	    final_df = pd.concat([final_df, df], axis=1, sort=True)
	final_df.to_csv(str(out_dir/'{}_rpkm_tables.csv'.format(run_name)), sep=',', index=True)
	np.log2(
		final_df+0.01).to_csv(out_dir/'log_tables.csv', 
		sep=',', 
		index=True,
		) 

def husky(log_file, out_dir):
	path_to_husky = str(
		Path(__file__).parent
		/'Algorithms'
		/'husky'
		/'TNBCtype_command_Latest.R'
		)
	# husky rscript would have to be changed to avoid changing dirs
	os.chdir(str(Path(__file__).parent/'Algorithms'/'husky'))
	sp.run([
		'R', 
		'--vanilla', 
		'-q', 
		'-f', 
		path_to_husky, 
		'--args', 
		log_file,
		out_dir,
		])

def lean(rpkm_table, out_dir):
	lean_path = str(
		Path(__file__).parent
		/'Algorithms'
		/'lean'
		/'rev_tnbc_lean.R')

	lean_dir = str(
		Path(__file__).parent
		/'Algorithms'
		/'lean')+'/'

	# Because lean adds .csv -- fix this at some point
	# Also change outfile name in Rscript
	rpkm_table = rpkm_table.split('.')[0]
	sp.run([
		'Rscript', 
		lean_path,
		lean_dir, 
		rpkm_table,
		])

def make_project(project_name):
	"""
	Makes analysis project.
	run_name = datetime.now().strftime('%H:%M:%S_%m-%d-%Y')
    when no argument for run_name is given?
    """

    # Analysis project structure
	master = Path(__file__).parent
	master_out = master / (project_name+'_analysis') 
	dirs = [str(master_out),											# project[0]
			str(master_out/'{}_bams'.format(project_name)),				# project[1]
			str(master_out/'{}_cufflinks'.format(project_name)),		# project[2]
			str(master_out/'{}_tables'.format(project_name))]			# project[3]]
	
	# Make sure project doesn't exit
	try:
		for d in dirs:
			Path(d).mkdir()
	except:
		a = str(input('Analysis project already exists, overwrite? (y/n)'))
		if a == 'y':
			for d in dirs:
				Path(d).mkdir(exist_ok=True)
			print('Analysis project successfully made!')
		else:
			print('Exiting...')
			sys.exit()
	# Returns list of dirs so paths can be used as i/o's for each function		
	return(dirs)

def qc_check(fastq_dir, out_dir):
	"""
	Produces a csv containing basecounts 
	and other relavent information about a fastq file.
	Can be used on one fastq or a directory of them.
	"""
	# Create list of fastq files, 
	if fastq_dir.endswith('.gz'):
		fastq_list = [Path(fastq_dir)]
	else: 
		fastq_list = [fastq for fastq in Path(fastq_dir).rglob('*.fastq.gz')]

	df = pd.DataFrame()

	# Append row to df for each fastq
	for fastq in fastq_list:
		try: 
	    	# Open compressed file
		    file = gzip.open(str(fastq), "rt")
		    # Make dictionary for counting 
		    d = collections.defaultdict(int)
		    # Read every 2nd line in fastq and count bases
		    for n, sequence in enumerate(file):
		        if n % 4 == 1:
		            for base in sequence:
		                d[base] += 1
		    file.close()
		    # Make OrderedDict to preserve desired order 
		    e = collections.OrderedDict([
		        ('filename', fastq.name),
		        ('total_read_count', d['\n']),
		        ('total_base_count', sum(d.values()) - d['\n']),
		        ('A', d['A']),
		        ('T', d['T']),
		        ('G', d['G']),
		        ('C', d['C']),
		        ('g_c_percentage', round((d['G']+d['C'])/(sum(d.values())-d['\n']), 2)),
		        ('N', d['N']),
		    ])
		    # Create a DataFrame from the OrderedDict and append to main DataFrame
		    temp_df = pd.DataFrame.from_dict(e, orient='index').T
		    df = pd.concat([df, temp_df], axis=0)
		    print('Finished parsing {}'.format(fastq.name))
		except:
			print('Problem with {}, moving on...'.format(fastq.name))
	# Export
	df.to_csv(str(Path(out_dir)/'qc.csv'), sep= ',', index=False)

def STAR(fastq_dir, out_dir):
	"""
	Runs STAR on fastq pairs contained in a given directory.
	- Can a posixpath object be used in sp without converting to string?
	- STAR doesn't allow for fastq path to contain commas, 
	it wont recognize the files, using os.chdir() to get around this.
	"""
	os.chdir(fastq_dir)
	dirs = Path(fastq_dir)
	run_name = dirs.name
	fastqs = sorted(list(dirs.glob('**/*fastq.gz')))

	fastq_list = [str(fastq).split(run_name)[1] 
					for fastq in fastqs 
					if not fastq.name.startswith('Undetermined')]
	
	n_pairs = int(len(fastq_list)/2)

	for i in range(0, n_pairs):
		R1 = fastq_list[2*i]
		R2 = fastq_list[2*i+1]
		bam_outfile = str(Path(out_dir)/Path(R1).name.split('.')[0])
		try:
			sp.run([
				star_path, 
				'--genomeDir', 
				index_path, 
				'--runThreadN', 
				n_proc, 
				'--readFilesCommand', 
				'zcat', 
				'--readFilesIn', 
				R1[1:], 
				R2[1:], 
				'--sjdbGTFfile', 
				gtf_path, 
				'--outFileNamePrefix', 
				bam_outfile, 
				'--outSAMtype', 
				'BAM', 
				'SortedByCoordinate'
				])		
		except:
			print("Error: {}".format(seq_name))
			continue 

def parse_adapters(fastq_1, fastq_2):
    name_1 = Path(fastq_1).name
    name_2 = Path(fastq_2).name
    print('Parsing {} and {}...'.format(name_1, Path(name_2)))
    fastqs = [fastq_1, fastq_2]
    d = collections.defaultdict(int)
    for fastq in fastqs:
        file = gzip.open(fastq, "rt")
        for n, sequence in enumerate(file):
            if n % 4 == 0:
                adapter = sequence.split(':')[-1][:-1]
                d[adapter] += 1
        file.close()
    print('Done parsing {} and {}...'.format(name_1, name_2))
    return(d)

def unidentified_qc(fastq_dir):
	# Assumes only 2 Unidentified fastqs in project dir
	adapter_csv = "/home/rick/Bioinf_Tools/adapterParse/RNAseq_index_list.csv"
	a = pd.read_csv(adapter_csv)
	d = Path(fastq_dir)
	c = collections.OrderedDict()
	fastqs = list(d.glob('**/*fastq.gz'))

	u = [str(fastq) for fastq in fastqs if str(fastq.name).startswith('Undetermined')]

	counts = parse_adapters(u[0], u[1])

	# Load c with IDs and sequences
	for i in range(15, 63):
	    adapter_ID = a.iloc[i,6]
	    adapter_sequence = a.iloc[i,5]
	    c[adapter_ID] = adapter_sequence

	final_dict = collections.OrderedDict([('Filename', str(Path(u[0]).name))])

	# Load adapter info, create dict (d) for all sequence/ID pairs
	for key in counts:
		adapter_IDs = [i for i, sequence in c.items() if sequence == key]
		for ID in adapter_IDs:
			final_dict[ID] = counts[key]

	final_dict['total_reads'] = sum(counts.values())
	temp_df = pd.DataFrame.from_dict(final_dict, orient='index')
	temp_df.to_csv(str(Path(u[0]).parent)+'/adapter_counts.csv', sep= ',', index=True)
	print('Done parsing {} and {}...'.format(u[0], u[1]))

def run_samples_through_ref_set(file):
	script = str(
		Path(__file__).parent
		/'ref_set_testing'
		/'run_samples_through_ref_set.py')

	sp.call([
		'python3',
		script,
		file])
	print('rpkm_merge_results with ref_set+lean is located in Bioinf_Tools/prnaSeq/ref_set_testing/tables')

if __name__ == '__main__':

	n_proc = str(multiprocessing.cpu_count())
	
	index_path = '/home/rick/Bioinf_Tools/STARIndex'
	star_path  = '/home/rick/STAR-2.5.2b/source/STAR'
	gtf_path   = '/home/rick/Bioinf_Tools/STARIndex/genes.gtf'

	parser = argparse.ArgumentParser()

	parser.add_argument('--all', nargs=1, help='Runs every process')
	parser.add_argument('--analyze', nargs=1, help='Performs analysis on fastq directory.')
	parser.add_argument('--qc', nargs=1, help='runs qc check on fastq dir, mostly used for testing')
	parser.add_argument('--quick', nargs=1, help = 'runs essentials without qc check for testing')
	args = vars(parser.parse_args())

	if args['all']:
		run_dir = Path(args['all'][0])
		run_name = run_dir.name
		project = make_project(run_name)

		# Main stuff
		bcl2fastq(args['all'][0])
		STAR(str(run_dir), project[1])
		cufflinks(project[1], project[2])
		cufflinks_merge(project[2], project[3], run_name)

		# Maybe check to see if log_tables is produced before moving on
		husky(str(Path(project[3])/'log_tables.csv'), project[3])
		#lean(str(Path(project[3])/'rpkm_tables.csv'.format(run_name)), project[3])

		# QC
		qc_check(str(run_dir), project[3])
		unidentified_qc(args['all'][0])
		run_samples_through_ref_set(project[3]+'/{}_rpkm_tables.csv'.format(run_name))

	if args['analyze']:
		run_dir = Path(args['analyze'][0])
		run_name = run_dir.name
		project = make_project(run_name)

		# Main stuff
		STAR(str(run_dir), project[1])
		cufflinks(project[1], project[2])
		cufflinks_merge(project[2], project[3], run_name)

		# Maybe check to see if log_tables is produced before moving on
		husky(str(Path(project[3])/'log_tables.csv'), project[3])
		#lean(str(Path(project[3])/'rpkm_tables.csv'.format(run_name)), project[3])

		# QC
		qc_check(str(run_dir), project[3])
		unidentified_qc(args['analyze'][0])
		run_samples_through_ref_set(project[3]+'/{}_rpkm_tables.csv'.format(run_name))
	if args['quick']:
		run_dir = Path(args['quick'][0])
		run_name = run_dir.name
		project = make_project(run_name)

		# Main stuff
		STAR(str(run_dir), project[1])
		cufflinks(project[1], project[2])
		cufflinks_merge(project[2], project[3], run_name)

		# Maybe check to see if log_tables is produced before moving on
		husky(str(Path(project[3])/'log_tables.csv'), project[3])
		#lean(str(Path(project[3])/'rpkm_tables.csv'.format(run_name)), project[3])

		# QC
		#qc_check(str(run_dir), project[3])
		#unidentified_qc(args['analyze'][0])
		run_samples_through_ref_set(project[3]+'/{}_rpkm_tables.csv'.format(run_name))	
	if args['qc']:
		run_dir = Path(args['qc'][0])
		run_name = run_dir.name
		project = make_project(run_name)
	pass
