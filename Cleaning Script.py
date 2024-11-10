# Import packages and assign aliases
import pandas as pd


# Read data and assign to variables
fpkm_df = pd.read_table('Data/fpkm1.txt', sep='\t')
gene_df = pd.read_table('Data/genes.tsv', sep='\t')
bud_df = pd.read_csv('Data/bud_counts.csv')


# Declare function to be used later
def dtype_strip(df, col_lst, dtype_lst):
    for col, dtype in zip(col_lst, dtype_lst):
        df[col] = df[col].astype(dtype)
        if dtype == 'string':
            df[col] = df[col].str.strip()
    return df


# fpkm_df Cleaning
fpkm_df = fpkm_df.melt(
    id_vars="time_points",
    var_name="time_point",
    value_name="fpkm1"
)
fpkm_df.rename(columns={'time_points': 'gene_transcript'}, inplace=True)
fpkm_cols = ['gene_transcript', 'time_point', 'fpkm1']
fpkm_dtypes = ['string', 'Int64', 'double']
fpkm_df = dtype_strip(fpkm_df, fpkm_cols, fpkm_dtypes)


# gene_df Cleaning
gene_cols = [
    'primary_identifier', 'secondary_identifier', 'symbol',
    'name', 'desc', 'qualifier',
    'chrom_identifier', 'chrom_start_pos',
    'chrom_end_pos', 'chrom_strand'
]
gene_df.columns = gene_cols
chrom_dict = {
    'chrI': 1, 'chrII': 2, 'chrIII': 3, 'chrV': 5, 'chrVIII': 8,
    'chrIX': 9, 'chrXI': 11, 'chrVI': 6, 'chrIV': 4, 'chrVII': 7,
    'chrX': 10, 'chrXII': 12, 'chrXIII': 13, 'chrXIV': 14, 'chrXV': 15,
    'chrXVI': 16, 'chrmt': 'mt', '<NA>': None
}
gene_df['chrom_identifier'] = gene_df['chrom_identifier'].replace(chrom_dict)

gene_dtypes = ['string', 'string', 'string', 'string', 'string', 'string',
               'string', 'Int64', 'Int64', 'Int64']
gene_df = dtype_strip(gene_df, gene_cols, gene_dtypes)


# bud_df Cleaning
bud_df = bud_df.iloc[0:51, 0:5]
bud_cols = ['time', 'zero_buds', 'one_bud', 'total_counted', 'percentage_budding']
bud_df.columns = bud_cols
bud_dtypes = ['Int64', 'Int64', 'Int64', 'Int64', 'double']
bud_df = dtype_strip(bud_df, bud_cols, bud_dtypes)


# Export Cleaned Dataframes
#gene_df.to_csv('genes_clean', index=False)
#fpkm_df.to_csv('fpkm1_clean', index=False)
#bud_df.to_csv('buds_clean', index=False)