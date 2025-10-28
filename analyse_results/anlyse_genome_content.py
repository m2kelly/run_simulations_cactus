import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
output_folder = '/home/maria/cactus_target_size/output_primates'
anlysis_file = '/home/maria/cactus_target_size/auxillary/primate_genome_heuristics2'
possible_species = ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']
'''
output_folder = '/home/maria/cactus_target_size/output_primates_extended'
analysis_file = '/home/maria/cactus_target_size/primate_extended/auxillaryprimate_genome_heuristics'
possible_species = ['Homo_sapiens','fullTreeAnc105','fullTreeAnc106','fullTreeAnc107','fullTreeAnc108','fullTreeAnc109','fullTreeAnc110']



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


speci_gc_percentage ={}
speci_non_syn_percentage = {}
speci_GC_ending_percentage = {}
for speci in possible_species:
    speci_df = pd.read_csv(f'{analysis_file}/{speci}', index_col =0)
    #restricting to cancer genes lists
    speci_df = speci_df[speci_df.index.isin(cancer_genes_list)]
    speci_sum_df = speci_df.sum()
    speci_gc_percentage[speci] = (speci_sum_df['C'] + speci_sum_df['G'])/(speci_sum_df['A'] + speci_sum_df['T']+speci_sum_df['C'] + speci_sum_df['G'])
    speci_non_syn_percentage[speci] = (speci_sum_df['M_n'])/(speci_sum_df['M_n'] + speci_sum_df['M_s'])
    speci_GC_ending_percentage[speci] = (speci_sum_df['syn_GC_ending'])/(speci_sum_df['syn_AT_ending'] + speci_sum_df['syn_GC_ending'])
    '''
    y = np.array([speci_sum_df[base] for base in bases])
    plt.figure()
    plt.pie(y, labels = bases)
    plt.title(speci)
    '''
    #plt.show() 

print(speci_gc_percentage)
point_sizes = {speci:30+i*30 for i,speci in enumerate(possible_species)}


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_gc_percentage[speci], label=speci, s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('gc content/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_content.svg', format='svg',bbox_inches='tight')

plt.figure()
plt.title('Cancer Genes')
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_non_syn_percentage[speci], label=speci, s = size)
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('M_n/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/c_M_n_content.svg', format='svg',bbox_inches='tight')


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_GC_ending_percentage[speci], label=speci, s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('GC_ending/total 4 fold degen sites')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_ending_4fold.svg', format='svg',bbox_inches='tight')


