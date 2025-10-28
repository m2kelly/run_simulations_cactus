import pandas as pd

#nonsyn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus/cactus/Anc4->hg38/nonsyn_output_data_preparation_cactus.txt'
#syn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus/cactus/Anc4->hg38/syn_output_data_preparation_cactus.txt'

nonsyn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_id/cactus/Anc4->hg38/nonsyn_output_data_preparation_cactus.txt'
syn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_id/cactus/Anc4->hg38/syn_output_data_preparation_cactus.txt'


cancer_genes = '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'

gene_id_df = pd.read_csv(gene_ids)

cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(gene_id_df, on=['gene_name'], how='inner')
cancer_genes_list = cancer_id_df['gene'].tolist()

# Read data
syn_data_df = pd.read_csv(syn_data_prep, index_col=0, delimiter='\t')
nonsyn_data_df = pd.read_csv(nonsyn_data_prep, index_col=0, delimiter='\t')
data_df = pd.merge(left = syn_data_df, right=nonsyn_data_df, how= 'inner', on = ['gene', 'l_m',	'l_k', 'l_s',	'm_obs',	'k_obs',	's_obs',	'L_gene','N_samples=1'])
data_df['muts'] = data_df['s_obs'] + data_df['m_obs'] + data_df['k_obs']

data_df['mu'] = data_df['lambda_s']/data_df['l_s']
data_df['l_n'] = data_df['l_m']+ data_df['l_k']
data_df['l_n_l_s'] = data_df['l_n']/ data_df['l_s']


cancer_df = data_df[data_df.index.isin(cancer_genes_list)]
non_cancer_df = data_df[~data_df.index.isin(cancer_genes_list)]

print(f'mean cancer mu {cancer_df['mu'].mean()}')
print(f'mean non cancer mu {non_cancer_df['mu'].mean()}')

print(f'mean cancer l_n {cancer_df['l_n'].mean()}')
print(f'mean non cancer l_n {non_cancer_df['l_n'].mean()}')

print(f'mean cancer l_n/l_s {cancer_df['l_n_l_s'].mean()}')
print(f'mean non cancer l_n/l_s {non_cancer_df['l_n_l_s'].mean()}')
cut_df = data_df[['lambda_n','lambda_s']]

#printing parameters
print(f'number of unique genes {len(cut_df.index.unique())}')
mut_df = data_df[data_df['muts']>0]
print(f'number of mutated genes {len(mut_df.index.unique())}')