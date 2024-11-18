# Import packages and assign aliases
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Read data, coerce dtypes, and assign to variables
integer_dict = {
    'chrom_start_pos' : 'Int64',
    'chrom_end_pos' : 'Int64',
    'chrom_strand' : 'Int64'
}
fpkm_df = pd.read_csv('Data/fpkm1_clean')
gene_df = pd.read_csv('Data/genes_clean', dtype=integer_dict)
bud_df = pd.read_csv('Data/buds_clean')

gene_df['name'] = gene_df['name'].astype('str')

# Hypothesis: 'The expression of mitochondrial replication genes shows less periodicity compared to the expression nuclear DNA replication genes'

# Define function to get gene data and assign complex names
def get_gene_complex(gene_df, gene_lst, complex_name):
    df = gene_df[gene_df['symbol'].isin(gene_lst)].copy()
    df.loc[:, 'complex'] = complex_name
    return df

# DNA Polymerases
pol_ep_lst = ['DPB2', 'DPB3', 'DPB4']
pol_ep = get_gene_complex(gene_df, pol_ep_lst, 'Pol Epsilon')

pol_al_lst = ['POL1', 'POL12', 'PRI2', 'PRI1']
pol_al = get_gene_complex(gene_df, pol_al_lst, 'Pol Alpha')

pol_de_lst = ['POL3', 'POL31', 'POL32']
pol_de = get_gene_complex(gene_df, pol_de_lst, 'Pol Delta')

