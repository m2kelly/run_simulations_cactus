from MutationTableGenerator import MutationTableGenerator

species = 'Anc4'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df_nocpg_nodup.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA'
N = 10000  #effective population size
Nb = 0.07
mutation_prob_file = '/home/maria/find_intron_matrix/output_cactus/rates'

generator = MutationTableGenerator(species, possible_species, exon_file, gene_annotations, output_dir,mutation_prob_file,N, Nb)
generator.run()

#then to load gene and seqs
'''
import pickle
with open("/home/maria/cactus_target_size/auxillary/non_syn_target_nocpg_anc/processed_gene_seqs.pkl", "rb") as f:
    while True:
        try:
            gene, seq = pickle.load(f)
            print(gene)
            print(seq)
            # process here
        except EOFError:
            break
'''

