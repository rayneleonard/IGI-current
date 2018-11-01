    #!/usr/bin/python3

# Runs on the directory containing R1/R2 output folders from CRC_MSI
# Grabs the msr_count.csv files recursively 
#
# $ python3 /path/to/script /path/to/input/dir /path/to/outfile.csv
# used pandas version 0.20.3

import pandas as pd
import glob
import os
import sys

# File list must contain r1/r2 R1/R2 in names and be composed of valid pairs
# Returns list of tuples containing r1/r2 pairs
def pair_by_end(files):
    R1 = [fq for fq in files if any(x in os.path.basename(fq) for x in ('r1', 'R1'))]
    R2 = [fq for fq in files if any(x in os.path.basename(fq) for x in ('r2', 'R2'))]
    if len(R1) == len(R2):
        return list(zip(sorted(R1),sorted(R2)))
    elif 0 in (len(R1), len(R2)):
        print('Empty fastq list')
    else:
        print('Problem with pairing')

if __name__ == '__main__':

    basedir = os.path.abspath(sys.argv[1])
    outfile = os.path.abspath(sys.argv[2])

    msrs = ['BAT25', 'BAT26', 'MONO27', 'NR21', 'NR24']
    files = glob.glob(basedir+'/*/*.csv')
    paired_files = pair_by_end(files)
    export_df = pd.DataFrame()

    for msr in msrs:
        df1 = pd.DataFrame()
        sums = []
        names = []

        for n, pair in enumerate(paired_files):
            R1, R2 = pair
            R1_name = os.path.basename(R1).strip('-msr_counts.csv')
            R2_name = os.path.basename(R2).strip('-msr_counts.csv')
            d1 = pd.DataFrame.from_csv(R1)[msr]
            d2 = pd.DataFrame.from_csv(R2)[msr]
            d3 = d1.add(d2)
            d3.index.name = "{} sum".format(R1_name)
            sums.append(d3)
            names.extend((R1_name, R2_name, d3.index.name))
            df1 = pd.concat([df1,d1,d2,d3], axis=1)

        df1.reset_index(inplace=True)
        df1.columns = range(df1.shape[1])
        max_df = pd.DataFrame(
            [['Sample ID', 'repeat size', 'count']+([None]*(df1.shape[1]-3))], 
            columns=range(df1.shape[1]), 
        )

        for i in sums:
            tdf = pd.DataFrame([[i.index.name, i.idxmax(), i.max()]])
            max_df = pd.concat([max_df, tdf])

        df1 = pd.concat([df1, max_df], ignore_index=True)
        df1.loc[-1] = [msr] + names
        df1 = df1.sort_index()
        df1.fillna('').head(30)
        export_df = pd.concat([export_df, df1])

    export_df.fillna('').to_csv(outfile, index=False, header=False)
