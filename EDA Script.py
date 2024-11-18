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

# Call function for nuclear DNA replication related genes
topoiso_lst = ['TOP1', 'TOP2', 'TOP3']
topiso = get_gene_complex(gene_df, topoiso_lst, 'Topoisomerase')

dna_heli_lst = ['HCS1', 'DNA2']
helicase = get_gene_complex(gene_df, dna_heli_lst, 'DNA Helicase')

ligase_lst = ['CDC9']
ligase = get_gene_complex(gene_df, ligase_lst, 'DNA Ligase')

pol_ep_lst = ['DPB2', 'DPB3', 'DPB4']
pol_ep = get_gene_complex(gene_df, pol_ep_lst, 'Pol Epsilon')

pol_al_lst = ['POL1', 'POL12', 'PRI2', 'PRI1']
pol_al = get_gene_complex(gene_df, pol_al_lst, 'Pol Alpha')

pol_de_lst = ['POL3', 'POL31', 'POL32']
pol_de = get_gene_complex(gene_df, pol_de_lst, 'Pol Delta')

telo_lst = ['EST1', 'EST2', 'EST3']
telo = get_gene_complex(gene_df, telo_lst, 'Telomerase')


# Concatenate nuclear DNA replication related gene DataFrames
nuclear_df = pd.concat([topiso, helicase, ligase, pol_ep, pol_al, pol_de, telo])



# Call function for mitochondrial replication related genes
# mt_ribo_lst = ['VAR1']
# mt_ribo = get_gene_complex(gene_df, mt_ribo_lst, 'Mitoribosomes')

pol_ga_lst = ['MIP1']
pol_ga = get_gene_complex(gene_df, pol_ga_lst, 'Pol Gamma')

mt_hel_lst = ['PIF1']
mt_hel = get_gene_complex(gene_df, mt_hel_lst, 'Mitochondrial Helicase')

# atp_syn_lst = [
#     "ATP1", "ATP3", "ATP7", "ATP16", "ATP5", "ATP17", "ATP12", 
#     "ATP2", "ATP14", "ATP10", "ATP11", "ATP4", "ATP15", "ATP20", 
#     "ATP18", "ATP8", "ATP6", "ATP19"
# ]
# atp_syn = get_gene_complex(gene_df, atp_syn_lst, 'ATP Synthase')

# cco_lst = [
#     "COX15", "COX6", "COX5B", "COX9", "COX20", "COX4", "COX13",
#     "COX18", "COX16", "COX17", "COX12", "COX8", "COX14", "COX7",
#     "COX5A", "COX11", "COX10", "COX19", "COX1", "COX2", "COX3"
# ]
# cco = get_gene_complex(gene_df, cco_lst, 'Cytochrome C Oxidase')

# Concatenate mitochondrial replication related gene DataFrames
mt_df = pd.concat([pol_ga, mt_hel])



# Assign gene values to a new column in each concatenated DataFrame
mt_df['genes'] =  'Mitochondria-Related'
nuclear_df['genes'] =  'Nuclear-Related'


# Concatenate both replication DataFrames
compar_df = pd.concat([nuclear_df, mt_df])

# Add time_point and fpkm1 values to each gene
test_df = pd.merge(
    compar_df, 
    fpkm_df, 
    how="inner", 
    left_on="symbol", 
    right_on="gene_transcript"
)



# Plot gene expression against time for each replication type
# plt.figure()

# fig = sns.lineplot(
#     x='time_point',
#     y=np.log10(test_df['fpkm1']),
#     hue='genes',
#     data=test_df
# )

# fig.set_title('Log10-Transformed Gene Expression vs Time')
# fig.set_xlabel('Time / mins')
# fig.set_ylabel('Log10(Gene Expression) / Log10(fpkm1)')
# fig.legend(title='Replication Genes')

# plt.show()



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
    legend=False
)

plt.suptitle('Log10-Transformed Gene Expression vs Time')
axes[0].set_xlabel('Time (mins)')
axes[1].set_xlabel('Time (mins)')
axes[0].set_ylabel('Log10(Gene Expression) (Log10(fpkm1))')
axes[0].set_title('nDNA Replication Genes')
axes[1].set_title('mtDNA Replication & Mitochondrial \n Component Genes')

plt.tight_layout()
plt.show()