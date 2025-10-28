from NeutralEvolutionSimulator import NeutralEvolutionSimulator
from CalculateSimulatedTarget import CalculateSimulatedTarget
from MultiplyTargets import MultiplyTargets
from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor

matrix_generator = MutationMatrixGenerator('', '', '', '', '')
print(f'calculating targets for simulation ')
target_caller = CalculateSimulatedTarget(
matrix_generator,
input_file = '/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/processed_gene_seqs.pkl',
output_dir = '/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/Anc4',
)
target_caller.run_parallel()


def process_one_signature(sig) -> str:
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/Anc4/*',
                            output_dir=f"/home/maria/run_simulations_cactus/output_anc4_hg38/Anc4",
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
with ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_one_signature, DEFAULT_SIGNATURE_LIST)





'''
#%%
#testing output
file = '/home/maria/run_simulations_cactus/output_anc4_hg38/simulated_genes_0.pkl'


import pandas as pd
import pickle
rows = []
with open(file, "rb") as f:
    while True:
        try:
            obj = pickle.load(f)
        except EOFError:
            break

        if isinstance(obj, list):     # came from chunked writer
            rows.extend(obj)
        else:                          # came from per-record writer
            rows.append(obj)
print(rows)
print(pd.DataFrame(rows,columns=['gene','seq','trinucs']))
'''
