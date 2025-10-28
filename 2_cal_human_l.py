from MutationMatrixGenerator import MutationMatrixGenerator
from MultiplyTargets import MultiplyTargets
from concurrent.futures import ProcessPoolExecutor

'''
species = 'hg38'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df_nocpg_nodup.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/hg38'


generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
generator.run()

'''
def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/hg38/*',
                            output_dir=f"/home/maria/run_simulations_cactus/output_anc4_hg38/hg38",
                            gene_strand_file='/home/maria/filter_transcripts/output/exon_merged_ids_strands')
    calc.calc_for_sim()
    return 

DEFAULT_SIGNATURE_LIST = [
    "SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d",
    "SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14",
    "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22a","SBS22b",
    "SBS23","SBS24","SBS25","SBS26","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
    "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40a","SBS40b","SBS40c","SBS41",
    "SBS42","SBS44","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91",
    "SBS92","SBS93","SBS94","SBS96","SBS97","SBS98","SBS99"
]

print('mulitplying by signatures')
with ProcessPoolExecutor(max_workers=10) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)

