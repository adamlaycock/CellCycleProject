# Import packages and assign aliases
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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

topoiso_lst = ['TOP1', 'TOP2', 'TOP3']
topoisomerases = gene_df[gene_df['symbol'].isin(topoiso_lst)]
topoisomerases['complex'] = 'Topoisomerase'

dna_heli_lst = ['HCS1']
helicase = gene_df[gene_df['symbol'].isin(dna_heli_lst)]
helicase['complex'] = 'DNA Helicase'

ligase_lst = ['DNL4']
ligase = gene_df[gene_df['symbol'].isin(ligase_lst)]
ligase['complex'] = 'DNA Ligase'

pol_ep_lst = ['DPB2', 'DPB3', 'DPB4']
pol_ep = gene_df[gene_df['symbol'].isin(pol_ep_lst)]
pol_ep['complex'] = 'Pol Epsilon'

pol_al_lst = ['POL1', 'POL12', 'PRI2', 'PRI1']
pol_al = gene_df[gene_df['symbol'].isin(pol_al_lst)]
pol_al['complex'] = 'Pol Alpha'

pol_de_lst = ['POL3', 'POL31', 'POL32']
pol_de = gene_df[gene_df['symbol'].isin(pol_de_lst)]
pol_de['complex'] = 'Pol Delta'

telo_lst = ['EST1', 'EST2', 'EST3']
telo = gene_df[gene_df['symbol'].isin(telo_lst)]
telo['complex'] = 'Telomerase'