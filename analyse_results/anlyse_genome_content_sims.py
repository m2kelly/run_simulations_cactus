import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
output_folder = '/home/maria/run_simulations_cactus/output_anc4_hg38'
analysis_file = '/home/maria/run_simulations_cactus/output_anc4_hg38'
possible_species = ['hg38',  'Anc4']
'''
#also chnage file reading line for sims vs actual 
output_folder = '/home/maria/run_simulations_cactus/output_anc4_hg38/genome_content'
analysis_file = '/home/maria/run_simulations_cactus/output_anc4_hg38'
possible_sims = [f'sim_{i}' for i in list(range(100))]

image_type='png'

plt.rcParams.update({'font.size': 18})

bases = ['A', 'C', 'G', 'T']

#extracting list of cancer gene ids 
cancer_genes = '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
gene_id_df = pd.read_csv(gene_ids)
cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(gene_id_df, on=['gene_name'], how='inner')
cancer_genes_list = cancer_id_df['gene'].tolist()



#anlaysing cancer gene 
speci_gc_percentage ={}
speci_non_syn_percentage = {}
speci_GC_ending_percentage = {}
for speci in possible_sims:
    #speci_df = pd.read_csv(f'{analysis_file}/{speci}_genome_content/{speci}', index_col =0)
    speci_df = pd.read_csv(f'{analysis_file}/{speci}_heuristics', index_col =0)
    
    
    #restricting to cancer genes lists
    speci_df = speci_df[speci_df.index.isin(cancer_genes_list)]
    speci_sum_df = speci_df.sum()
    speci_gc_percentage[speci] = (speci_sum_df['C'] + speci_sum_df['G'])/(speci_sum_df['A'] + speci_sum_df['T']+speci_sum_df['C'] + speci_sum_df['G'])
    speci_non_syn_percentage[speci] = (speci_sum_df['M_n'])/(speci_sum_df['M_n'] + speci_sum_df['M_s'])
    speci_GC_ending_percentage[speci] = (speci_sum_df['syn_GC_ending'])/(speci_sum_df['syn_AT_ending'] + speci_sum_df['syn_GC_ending'])
  

print(speci_gc_percentage) 
print(np.mean(list(speci_gc_percentage.values()))) 

'''
#analysing non cancer genes
speci_gc_percentage ={}
speci_non_syn_percentage = {}
speci_GC_ending_percentage = {}
for speci in possible_species:
    speci_df = pd.read_csv(f'{analysis_file}/{speci}_genome_content/{speci}', index_col =0)
    #restricting to cancer genes lists
    speci_df = speci_df[~speci_df.index.isin(cancer_genes_list)]
    speci_sum_df = speci_df.sum()
    speci_gc_percentage[speci] = (speci_sum_df['C'] + speci_sum_df['G'])/(speci_sum_df['A'] + speci_sum_df['T']+speci_sum_df['C'] + speci_sum_df['G'])
    speci_non_syn_percentage[speci] = (speci_sum_df['M_n'])/(speci_sum_df['M_n'] + speci_sum_df['M_s'])
    speci_GC_ending_percentage[speci] = (speci_sum_df['syn_GC_ending'])/(speci_sum_df['syn_AT_ending'] + speci_sum_df['syn_GC_ending'])
    
'''

