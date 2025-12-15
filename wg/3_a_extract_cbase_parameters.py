import pandas as pd
output_file = '/home/maria/run_simulations_cactus/auxiliary_wg/cbase.csv'
#nonsyn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus/cactus/Anc4->hg38/nonsyn_output_data_preparation_cactus.txt'
#syn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus/cactus/Anc4->hg38/syn_output_data_preparation_cactus.txt'

nonsyn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_nocpg/cactus/Anc4->hg38/nonsyn_output_data_preparation_cactus.txt'
syn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_nocpg/cactus/Anc4->hg38/syn_output_data_preparation_cactus.txt'

#nonsyn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_id/cactus/Anc4->hg38/nonsyn_output_data_preparation_cactus.txt'
#syn_data_prep = '/home/maria/hossam_cbase_exons/CBase_external_matrix/Output_cactus_id/cactus/Anc4->hg38/syn_output_data_preparation_cactus.txt'


# Read data
syn_data_df = pd.read_csv(syn_data_prep, index_col=0, delimiter='\t')
nonsyn_data_df = pd.read_csv(nonsyn_data_prep, index_col=0, delimiter='\t')
data_df = pd.merge(left = syn_data_df, right=nonsyn_data_df, how= 'inner', on = ['gene', 'l_m',	'l_k', 'l_s',	'm_obs',	'k_obs',	's_obs',	'L_gene','N_samples=1'])
data_df['muts'] = data_df['s_obs'] + data_df['m_obs'] + data_df['k_obs']

cut_df = data_df[['lambda_n','lambda_s']]
cut_df.to_csv(output_file)

#printing parameters
print(f'number of mutations {data_df['muts'].sum()}')
print(f'number of unique genes {len(cut_df.index.unique())}')
mut_df = data_df[data_df['muts']>0]
print(f'number of mutated genes {len(mut_df.index.unique())}')