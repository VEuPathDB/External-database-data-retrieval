from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import pandas as pd
import glob

def get_data(fasta):
    fasta_data = [['file_name','contig','length','checksum']]
    file_name = fasta
    for record in SeqIO.parse(fasta,'fasta'):
        contig = record.id
        length = len(record.seq)
        checksum = seguid(record.seq.upper())
        fasta_data.append([file_name,contig,length,checksum])
    df = pd.DataFrame(fasta_data[1:], columns=fasta_data[0])
    return df

# find fasta files and run    
files = glob.glob(r"*.fasta", recursive=True)
all_data_df = pd.DataFrame(columns=['file_name','contig','length','checksum'])
for file in files:
    fasta_df = get_data(file)
    all_data_df = all_data_df.merge(fasta_df, how='outer')

# tidying
all_data_df['database'] = 0
all_data_df.loc[all_data_df['file_name'].str.contains('FungiDB'),'database'] = 'veupath'
all_data_df.loc[~all_data_df['file_name'].str.contains('FungiDB'),'database'] = 'CGD'
all_data_df.loc[all_data_df['file_name'].str.contains('POMBASE'),'database'] = 'PomBase'
all_data_df.loc[all_data_df['file_name'].str.contains('SGD'),'database'] = 'SGD'
all_data_df.loc[all_data_df['file_name'].str.contains('GCA_'),'database'] = 'GenBank'
all_data_df.loc[all_data_df['file_name'].str.contains('GCF_'),'database'] = 'RefSeq'
all_data_df['organism'] = all_data_df['file_name'].str.replace('_chromosomes.fasta','')\
    .str.replace('_current','').str.replace('_Genome.fasta','').str.replace('FungiDB-65_','')\
        .str.replace('_','').str.replace('-','')
all_data_df.sort_values(by=['organism','length','checksum'])