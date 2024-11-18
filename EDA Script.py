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

# Nuclear DNA Replicated Related Genes
topoiso_lst = [
    'TOP1', 
    'TOP2', 
    'TOP3'
]
topiso = gene_df[gene_df['symbol'].isin(topoiso_lst)].copy()
topiso.loc[:, 'complex'] = 'Topoisomerase'

dna_heli_lst = [
    'HCS1'
]
helicase = gene_df[gene_df['symbol'].isin(dna_heli_lst)].copy()
helicase.loc[:, 'complex'] = 'DNA Helicase'

ligase_lst = [
    'CDC9'
]
ligase = gene_df[gene_df['symbol'].isin(ligase_lst)].copy()
ligase.loc[:, 'complex'] = 'DNA Ligase'

pol_ep_lst = [
    'DPB2', 
    'DPB3', 
    'DPB4'
]
pol_ep = gene_df[gene_df['symbol'].isin(pol_ep_lst)].copy()
pol_ep.loc[:, 'complex'] = 'Pol Epsilon'

pol_al_lst = [
    'POL1', 
    'POL12', 
    'PRI2', 
    'PRI1'
]
pol_al = gene_df[gene_df['symbol'].isin(pol_al_lst)].copy()
pol_al.loc[:, 'complex'] = 'Pol Alpha'

pol_de_lst = [
    'POL3', 
    'POL31', 
    'POL32'
]
pol_de = gene_df[gene_df['symbol'].isin(pol_de_lst)].copy()
pol_de.loc[:, 'complex'] = 'Pol Delta'

telo_lst = [
    'EST1', 
    'EST2', 
    'EST3'
]
telo = gene_df[gene_df['symbol'].isin(telo_lst)].copy()
telo.loc[:, 'complex'] = 'Telomerase'

nuclear_df = pd.concat([topiso, helicase, ligase, pol_ep, pol_al, pol_de, telo])

# Mitochondrial Replication Related Genes
mt_ribo_lst = [
    'VAR1'
]
mt_ribo = gene_df[gene_df['symbol'].isin(mt_ribo_lst)].copy()
mt_ribo.loc[:, 'complex'] = 'Mitoribosomes'

pol_ga_lst = [
    'MIP1'
]
pol_ga = gene_df[gene_df['symbol'].isin(pol_ga_lst)].copy()
pol_ga.loc[:, 'complex'] = 'Pol Gamma'

atp_syn_lst = [
    "ATP1", "ATP3", "ATP7", "ATP16", "ATP5", "ATP17", "ATP12", 
    "ATP2", "ATP14", "ATP10", "ATP11", "ATP4", "ATP15", "ATP20", 
    "ATP18", "ATP8", "ATP6", "ATP19"
]
atp_syn = gene_df[gene_df['symbol'].isin(atp_syn_lst)].copy()
atp_syn.loc[:, 'complex'] = 'ATP Synthase'

cco_lst = [
    "COX15", "COX6", "COX5B", "COX9", "COX20", "COX4", "COX13",
    "COX18", "COX16", "COX17", "COX12", "COX8", "COX14", "COX7",
    "COX5A", "COX11", "COX10", "COX19", "COX1", "COX2", "COX3"
]
cco = gene_df[gene_df['symbol'].isin(cco_lst)].copy()
cco.loc[:, 'complex'] = 'Cytochrome C Oxidase'

mt_df = pd.concat([mt_ribo, pol_ga, atp_syn, cco])

mt_df['genes'] =  'Mitochondria-Related'
nuclear_df['genes'] =  'Nuclear-Related'

compar_df = pd.concat([nuclear_df, mt_df])

test_df = pd.merge(compar_df, fpkm_df, how="inner", left_on="symbol", right_on="gene_transcript")

plt.figure()

sns.lineplot(
    x='time_point',
    y=np.log(test_df['fpkm1']),
    hue='genes',
    data=test_df
)

plt.show()