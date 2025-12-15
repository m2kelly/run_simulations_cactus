from MutationMatrixGenerator import MutationMatrixGenerator


species = 'Anc4'
possible_species = ['hg38', 'Anc4']
exon_file = '/home/maria/run_simulations_cactus/auxillary_wg/anc4_nodup_nocpg.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/Anc4_wg'


generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
generator.run()