import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


output_folder = '/home/maria/run_simulations_cactus/output_anc4_hg38'
analysis_file = '/home/maria/run_simulations_cactus/output_anc4_hg38'
possible_species = ['hg38',  'Anc4']
'''
#also chnage file reading line for sims vs actual 
output_folder = '/home/maria/run_simulations_cactus/output_anc4_hg38/genome_content'
analysis_file = '/home/maria/run_simulations_cactus/output_anc4_hg38'
possible_species = ['sim_0']
'''
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
for speci in possible_species:
    speci_df = pd.read_csv(f'{analysis_file}/{speci}_genome_content/{speci}', index_col =0)
    #speci_df = pd.read_csv(f'{analysis_file}/{speci}_heuristics', index_col =0)
    
    print(speci_df.columns)
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
print(f'cancer gc Delta from anc to human {(speci_gc_percentage['Anc4']-speci_gc_percentage['hg38'])/(speci_gc_percentage['Anc4']+speci_gc_percentage['hg38'])}')

plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_gc_percentage[speci], label=speci, s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('gc content/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_content.{image_type}', format=image_type,bbox_inches='tight')

plt.figure()
plt.title('Cancer Genes')
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_non_syn_percentage[speci], label=speci, s = size)
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('M_n/valid aligned bases')
ax = plt.gca()
# 💡 This disables both scientific notation and the "+7.746e" style offset text
ax.ticklabel_format(style='plain', axis='y', useOffset=False)
ax.yaxis.get_major_formatter().set_scientific(False)
ax.yaxis.get_major_formatter().set_useOffset(False)
#plt.legend()
plt.savefig(f'{output_folder}/c_M_n_content.{image_type}', format=image_type,bbox_inches='tight')


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_GC_ending_percentage[speci], label=speci, s = size)
plt.title('Cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('GC_ending/total 4 fold degen sites')
#plt.legend()
plt.savefig(f'{output_folder}/c_gc_ending_4fold.{image_type}', format=image_type,bbox_inches='tight')




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
    y = np.array([speci_sum_df[base] for base in bases])
    plt.figure()
    plt.pie(y, labels = bases)
    plt.title(speci)
    '''
    #plt.show() 

print(speci_gc_percentage)
point_sizes = {speci:30+i*30 for i,speci in enumerate(possible_species)}
print(f'non cancer gc Delta from anc to human {(speci_gc_percentage['Anc4']-speci_gc_percentage['hg38'])/(speci_gc_percentage['Anc4']+speci_gc_percentage['hg38'])}')


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_gc_percentage[speci], label=speci, s = size)
plt.title('Non cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('gc content/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/nc_gc_content.{image_type}', format=image_type,bbox_inches='tight')

print(speci_non_syn_percentage)
plt.figure()
plt.title('Non cancer Genes')
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_non_syn_percentage[speci], label=speci, s = size)
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('M_n/valid aligned bases')
#plt.legend()
plt.savefig(f'{output_folder}/nc_M_n_content.{image_type}', format=image_type,bbox_inches='tight')


plt.figure()
for speci in possible_species:
    size = point_sizes.get(speci, 40)
    plt.scatter(speci, speci_GC_ending_percentage[speci], label=speci, s = size)
plt.title('Non cancer Genes')
plt.xlabel('species')
plt.xticks(rotation=70)
plt.ylabel('GC_ending/total 4 fold degen sites')
#plt.legend()
plt.savefig(f'{output_folder}/nc_gc_ending_4fold.{image_type}', format=image_type,bbox_inches='tight')

