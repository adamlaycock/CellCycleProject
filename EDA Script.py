# Import packages and assign aliases
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

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

# Hypothesis: 'The average expression of mtDNA polymerase genes is lower than the average expression of nDNA polymerase genes'
# Null Hypothesis: 'There is no difference in average expression of mtDNA polymerase genes compared to nDNA polymerase genes'

# Define function to get gene data and assign complex names
def get_gene_complex(gene_df, gene_lst, complex_name):
    df = gene_df[gene_df['symbol'].isin(gene_lst)].copy()
    df.loc[:, 'complex'] = complex_name
    return df

# Nuclear DNA Polymerases
pol_ep_lst = ['DPB2', 'DPB3', 'DPB4']
pol_ep = get_gene_complex(gene_df, pol_ep_lst, 'Pol Epsilon')

pol_al_lst = ['POL1', 'POL12', 'PRI2', 'PRI1']
pol_al = get_gene_complex(gene_df, pol_al_lst, 'Pol Alpha')

pol_de_lst = ['POL3', 'POL31', 'POL32']
pol_de = get_gene_complex(gene_df, pol_de_lst, 'Pol Delta')


# Mitochondrial DNA Polymerases
pol_ga_lst = ['MIP1']
pol_ga = get_gene_complex(gene_df, pol_ga_lst, 'Pol Gamma')

mtDNA_df = pol_ga
nDNA_df = pd.concat([pol_ep, pol_al, pol_de])


# Assign gene values to a new column in each concatenated DataFrame
mtDNA_df['genes'] =  'Mitochondria-Related'
nDNA_df['genes'] =  'Nuclear-Related'


# Concatenate both replication DataFrames
compar_df = pd.concat([nDNA_df, mtDNA_df])

# Add time_point and fpkm1 values to each gene
test_df = pd.merge(
    compar_df, 
    fpkm_df, 
    how="inner", 
    left_on="symbol", 
    right_on="gene_transcript"
)

plt.figure()

fig, axes = plt.subplots(1, 2, sharey=True)

nuc_df = test_df[test_df['genes'] == 'Nuclear-Related']
mta_df = test_df[test_df['genes'] == 'Mitochondria-Related']

fig = sns.lineplot(
    x='time_point',
    y=np.log10(test_df['fpkm1']),
    data=nuc_df,
    ax=axes[0],
    legend=False
)

fig = sns.lineplot(
    x='time_point',
    y=np.log10(test_df['fpkm1']),
    data=mta_df,
    ax=axes[1],
    legend=False,
    color='orange'
)

plt.suptitle('Log10-Transformed Gene Expression vs Time')
axes[0].set_xlabel('Time (mins)')
axes[1].set_xlabel('Time (mins)')
axes[0].set_ylabel('Log10(Gene Expression) (Log10(fpkm1))')
axes[0].set_title('nDNA Polymerase Genes')
axes[1].set_title('mtDNA Polymerase Genes')

plt.tight_layout()
plt.show()

nDNA_fpkm1 = nuc_df['fpkm1']
mtDNA_fpkm1 = mta_df['fpkm1']

U, p = mannwhitneyu(nDNA_fpkm1, mtDNA_fpkm1)

# Average gene expression in mtDNA and nDNA were 12.59FPKm1 and 34.24FPKm1; 
# the distributions showed a significant difference 
# (Mann–Whitney U=22369.5, n1=50, n2=500, p=3.22e-20).