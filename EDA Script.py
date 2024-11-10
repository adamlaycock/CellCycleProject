# Import packages and assign aliases
import pandas as pd

# Read data, coerce dtypes, and assign to variables
integer_dict = {
    'chrom_start_pos' : 'Int64',
    'chrom_end_pos' : 'Int64',
    'chrom_strand' : 'Int64'
}
fpkm_df = pd.read_csv('Data/fpkm1_clean')
gene_df = pd.read_csv('Data/genes_clean', dtype=integer_dict)
bud_df = pd.read_csv('Data/buds_clean')