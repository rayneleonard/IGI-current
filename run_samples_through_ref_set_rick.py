'''run sample + refset through lean'''
## python3 /home/summer/Desktop/ref_set_testing/run_samples_through_ref_set.py ~/Desktop/rpkm.csv
## don't fucking shortcut the above path. WILL BREAK. 
## clear out storage folder

from prnaSeq import lean
from pathlib import Path

import pandas as pd
import os
import sys
import glob
import shutil

PROJECT_ROOT = os.path.dirname(__file__)

outdir   = PROJECT_ROOT+"/tables/"
lean_out = PROJECT_ROOT+"/lean_output/"
storage  = PROJECT_ROOT+"/ref_set_with_samples/"
refset   = PROJECT_ROOT+"/master_ref_set.csv"

# 102 Genes
genes = ['UBD','CXCL13','CXCL11','IDO1','IGHM','CCL5',
        'CXCL10','NTN3','GZMB','IGKC','LOC652493',
        'PSMB9','GBP1','FGL2','IL23A','CD52',
        'RARRES3','IGJ','RTP4','TNFSF10','ITM2A',
        'TNFAIP8','GABRP','KYNU','S100A8','PTGDS',
        'PMAIP1','FABP7','KMO','XBP1','CHI3L1',
        'DEFB1','DUSP5','ART3','NEK2','C7',
        'TSC22D3','CDC20','SPTLC2','MYBL1','IL33',
        'SIDT1','LBP','GALNT7','JAM2','S100A7',
        'PI3','GCHFR','PLEKHB1','MFAP4','TRIM68',
        'TCF7L1','INPP4B','FOXA1','SPRR1B','FGFR4',
        'ADRA2A','SFRP1','ABCA8','FASN','CRAT',
        'KCNK5','ALDH3B2','MID1','CYP4F8','DHCR24',
        'IGFBP4','HGD','LHFP','BLVRB','SPDEF',
        'CROT','GPR87','ZCCHC24','THBS4','SEMA3C',
        'OGN','KRT6A','COL5A2','DBI','AKAP12','HTRA1',
        'FOXC1','CD36','MIA','S100A1','TFAP2B',
        'KIAA1324','SERHL2','COL5A1','UGT2B28','KRT16',
        'KRT6B','SRPX','SOX10','AZGP1','APOD',
        'KRT17','KRT14','COL2A1','ASPN','SCRG1']

# Read in count table
rpkm = sys.argv[1] 
df = pd.read_csv(rpkm)

# Filter down to 102 genes
filtered_df = df[df['gene_short_name'].isin(genes)]

# Make df containing samples to be pushed through refset
sample_df = filtered_df.reset_index()
sample_df = sample_df[sample_df.columns[2:]]

# Make refset and index df
rdf = pd.read_csv(refset)
idf = rdf['gene_short_name']
rdf = rdf[rdf.columns[2:]]

csv_list = []

for sample in sample_df:
    print(sample)
    final_df = pd.concat([idf, sample_df[sample], rdf], axis=1)
    final_df.to_csv(outdir+"{}_RS.csv".format(sample), index=False)
    csv_list.append(outdir+"{}_RS.csv".format(sample))

for csv in csv_list:
	lean(csv, outdir)
	shutil.move(csv, storage+str(Path(csv).name))

# Read in csv list from directory
csv_list = glob.glob(outdir+'*log.csv')

# Creating DataFrame to append columns to
final_df = pd.DataFrame()

for csv in csv_list:
    df = pd.read_csv(csv, index_col=0).transpose()
    col_0 = df[df.columns[0]]
    # Formatting magic, i.e. "V515.06187.A1" -> "V515.06187.A1_results_centroids_(n)"
    file_name = csv.split("/")[-1]
    result_number = file_name.split('.')[0]
    result_number = '_'.join(result_number.split('_')[-2:])
    sample_name = '_'.join([df.columns[1], result_number])
    df.rename(columns={df.columns[1]: sample_name}, inplace=True)
    # Concat each sample to final dataframe
    final_df = pd.concat([final_df, df[df.columns[0]]], axis=1)
    shutil.move(csv, lean_out+str(Path(csv).name))
    
final_df = pd.concat([col_0,final_df], axis=1).T
final_df = final_df.fillna('NA')

# This puts columns in the correct order
cols = ['Base_Call', 'Strong_Call', 'IM_Call_Vandy', 
        'Confounding_Class', 'BL1_corr', 'BL2_corr', 
        'LAR_corr', 'M_corr', 'MSL_corr', 'IM_corr', 
        'BL1_pvalue', 'BL2_pvalue', 'LAR_pvalue', 
        'M_pvalue', 'MSL_pvalue', 'IM_pvalue']

final_df = final_df[cols]

outname = str(Path(rpkm).name).split('.')[0]
final_df.to_csv(outdir+outname+'_merge_results.csv')
