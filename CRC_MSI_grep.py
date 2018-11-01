import matplotlib
matplotlib.use('Agg')
import os
import sys
import glob
import filecmp
import subprocess as sp
import pandas as pd
import matplotlib.pyplot as plt

project_dir = os.path.dirname(os.path.abspath(__file__))

MSRS = {
    'BAT25': { 
        'LFR': {
            'PATTERN': 'GAGTTTTGTGTTTTGTTTTTTTGA', 
            'NUC': 'T',
            'SERIES': ['A{}G'.format('T'*i) for i in range(10,31)],
        },
        'RFR': {
            'PATTERN': 'CTAAAATGCTCTGTTCTC', 
            'NUC': 'A',
            'SERIES': ['C{}T'.format('A'*i) for i in range(10,31)],
        },
    },
    'BAT26': { 
        'LFR': {
            'PATTERN': 'CTTAACCTTTTTCAGGT', 
            'NUC': 'A',
            'SERIES': ['T{}G'.format('A'*i) for i in range(10,31)],
        },
        'RFR': {
            'PATTERN': 'TCAACATTTTTAACCC', 
            'NUC': 'T',
            'SERIES': ['C{}A'.format('T'*i) for i in range(10,31)],
        },
    },
    'NR21': {
        'LFR': {
            'PATTERN': 'AAAAGTGTTGCT', 
            'NUC': 'A',
            'SERIES': ['T{}G'.format('A'*i) for i in range(10,31)],
        },
        'RFR': {
            'PATTERN': 'ATGTATGTCTCCCCTGGCC', 
            'NUC': 'T',
            'SERIES': ['C{}A'.format('T'*i) for i in range(10,31)],
        },
    },
    'NR24': {
        'LFR': {
            'PATTERN': 'CCCAGTCCTA', 
            'NUC': 'T',
            'SERIES': ['A{}G'.format('T'*i) for i in range(10,31)],
        },
        'RFR': {
            'PATTERN': 'AGACTCTGTCTCAC', 
            'NUC': 'A',
            'SERIES': ['C{}T'.format('A'*i) for i in range(10,31)],
        },
    },
    'MONO27': {
        'LFR': {
            'PATTERN': 'GTAAAACCAGGA', 
            'NUC': 'T',
            'SERIES': ['A{}G'.format('T'*i) for i in range(10,31)],
        },
        'RFR': {
            'PATTERN': 'AGCAAGACTCTGCCTC', 
            'NUC': 'A',
            'SERIES': ['C{}C'.format('A'*i) for i in range(10,31)],
        },
    },
}

def get_fastq_name(fastq):
    name = fastq.split('/')[-1]
    return name.split('.')[0]
    
def is_same_file(file_one, file_two):
    return filecmp.cmp(file_one, file_two)

def seqtk(fastq, outfile):
    arg = "seqtk seq -a {} > {}".format(fastq, outfile)
    sp.check_call(arg, shell=True)

def grep(file, pattern, out_file, c=False):
    with open(out_file, 'w') as of:
        with open(file, 'r') as f:
            if c == True:
                count = 0
                for line in f:
                    if pattern in line:
                        count += 1
                of.write(str(count))
            else:
                for line in f:
                    if pattern in line:
                        of.write(line)

def agrep(file, pattern, out_file):
    arg = "grep -c '{}' {} > {}".format(pattern, file, out_file)
    sp.check_call(arg, shell=True)

def make_seq_dirs(fastq):
    fastq_name = get_fastq_name(fastq)
    project = {}
    d = {msr+name: project_dir+'/{}/{}'.format(fastq_name, msr+name) 
         for msr in MSRS 
         for name in ['-out', '-repeat']}
    paths = [v for k,v in d.items()]
    for p in paths:
        os.makedirs(p, exist_ok=True)
    d['ROOT'] = project_dir + '/{}'.format(fastq_name)
    project[fastq_name] = d
    return project

def get_file_contents(file):
    with open(file, 'r') as f:
        return int(f.read())
    
def analyze_fastq_pair(R1, R2):
    final_dict = {}
    root_paths = []

    for i in [R1, R2]:
        p = make_seq_dirs(i)
        name = get_fastq_name(i)

        root = p[name]['ROOT']
        root_paths.append(root)

        fasta = p[name]['ROOT'] + '/{}.fasta'.format(name)
        A10 = p[name]['ROOT'] + '/{}-A10.txt'.format(name)
        T10 = p[name]['ROOT'] + '/{}-T10.txt'.format(name)
        seqtk(i, fasta)
        grep(fasta, 'AAAAAAAAAA', A10)
        grep(fasta, 'TTTTTTTTTT', T10)
        
        second_temp_dict = {}

        for msr, d in MSRS.items():
            A_out = p[name][msr+'-out'] + '/{}-{}.txt'.format(name+'-A10', msr)
            T_out = p[name][msr+'-out'] + '/{}-{}.txt'.format(name+'-T10', msr)
            temp_dict = {}

            if i == R1:
                if d['LFR']['NUC'] == 'A':
                    grep(A10, d['LFR']['PATTERN'], A_out)
                    for n, pattern in enumerate(d['LFR']['SERIES']):
                        out_file = p[name][msr+'-repeat']+'/{}-{}-{}-{}-repeats.txt'.format(name, msr, n+10, 'A')
                        try:
                            agrep(A_out, pattern, out_file)
                            n_counts = get_file_contents(out_file)
                            temp_dict[n+10] = n_counts
                        except sp.CalledProcessError:
                            pass

                elif d['LFR']['NUC'] == 'T':
                    grep(T10, d['LFR']['PATTERN'], T_out)
                    for n, pattern in enumerate(d['LFR']['SERIES']):
                        out_file = p[name][msr+'-repeat']+'/{}-{}-{}-{}-repeats.txt'.format(name, msr, n+10, 'T')
                        try:
                            agrep(T_out, pattern, out_file)
                            n_counts = get_file_contents(out_file)
                            temp_dict[n+10] = n_counts
                        except sp.CalledProcessError:
                            pass

            elif i == R2:
                if d['RFR']['NUC'] == 'A':
                    grep(A10, d['RFR']['PATTERN'], A_out)
                    for n, pattern in enumerate(d['RFR']['SERIES']):
                        out_file = p[name][msr+'-repeat']+'/{}-{}-{}-{}-repeats.txt'.format(name, msr, n+10, 'A')
                        try:
                            agrep(A_out, pattern, out_file)
                            n_counts = get_file_contents(out_file)
                            temp_dict[n+10] = n_counts
                        except sp.CalledProcessError:
                            pass           

                if d['RFR']['NUC'] == 'T':
                    grep(T10, d['RFR']['PATTERN'], T_out)
                    for n, pattern in enumerate(d['RFR']['SERIES']):
                        out_file = p[name][msr+'-repeat']+'/{}-{}-{}-{}-repeats.txt'.format(name, msr, n+10, 'T')
                        try:
                            agrep(T_out, pattern, out_file)
                            n_counts = get_file_contents(out_file)
                            temp_dict[n+10] = n_counts
                        except sp.CalledProcessError:
                            pass

            second_temp_dict[msr] = temp_dict            
        final_dict[i] = second_temp_dict   
    return final_dict, root_paths

def create_bar_plot_and_csv(fastq, count_dict, out_dir):
        name = get_fastq_name(fastq)
        df = pd.DataFrame.from_dict(count_dict[fastq]).fillna(0)
        #df.index.name = os.path.basename(fastq)
        df.to_csv(out_dir+'/{}-msr_counts.csv'.format(name))
        plt.figure()
        plt.rcParams["figure.figsize"] = (15,8)
        # stacked=True parameter gave weird numbers
        a = df.plot(kind='bar', alpha=0.75)
        plt.xlabel('Repeats', fontsize=12)
        plt.ylabel('Counts', fontsize=12)
        plt.title('{}'.format(name), fontsize=14, fontweight='bold')
        plt.xticks(rotation=0)
        plt.savefig(out_dir+'/{}.png'.format(name))
                            
if __name__ == '__main__':
    fastq_dir = os.path.abspath(sys.argv[1])
    fastqs = sorted(glob.glob(fastq_dir+'/*fastq.gz'))
    n_pairs = int(len(fastqs)/2)

    for i in range(0, n_pairs):
        R1 = fastqs[2*i]
        R2 = fastqs[2*i+1]
        
        count_dict, root_paths = analyze_fastq_pair(R1=R1, R2=R2)

        # create tuple of fastq and its project path
        t = list(zip([R1, R2], root_paths))

        for fastq, base_path in t:
            create_bar_plot_and_csv(fastq, count_dict, base_path)

