from NeutralEvolutionSimulator import NeutralEvolutionSimulator
from CalculateSimulatedTarget import CalculateSimulatedTarget
from MultiplyTargets import MultiplyTargets
from MutationMatrixGenerator import MutationMatrixGenerator
from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import shutil

sim_start_no = 30
no_of_sims = 100

def process_one_signature(args) -> str:
    sig,x = args
    calc = MultiplyTargets(signature=sig,
                            input_glob= f'/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/sim_{x}/*',
                            output_dir=f"/home/maria/run_simulations_cactus/output_anc4_hg38/sim_{x}",
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
matrix_generator = MutationMatrixGenerator('', '', '', '', '')

for x in range(sim_start_no,no_of_sims):
    print(f'running simulation {x}')
    simulator = NeutralEvolutionSimulator(
        n_subs=72047,
        start_speci_file='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/processed_gene_seqs.pkl',
        cbase_output='/home/maria/run_simulations_cactus/auxiliary/cbase.csv',
        syn_probs_dir='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/syn_target',
        non_syn_probs_dir='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/nonsyn_target',
        output=f'/home/maria/run_simulations_cactus/output_anc4_hg38/simulated_genes_{x}.pkl',
        gene_name_map = '/home/maria/filter_transcripts/output/gene_name_id.csv')
    simulator.run()
    print(f'calculating targets for simulation {x}')
    target_caller = CalculateSimulatedTarget(
    matrix_generator,
    input_file = f'/home/maria/run_simulations_cactus/output_anc4_hg38/simulated_genes_{x}.pkl',
    output_dir = f'/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/sim_{x}',
    )
    target_caller.run_parallel()
    print(f'mutiplying targets for simulation {x}')
    args_list = [(sig, x) for sig in DEFAULT_SIGNATURE_LIST]
    with ProcessPoolExecutor(max_workers=15) as executor:
        executor.map(process_one_signature, args_list)
    p=Path(f'/home/maria/run_simulations_cactus/target_sizes_anc4_hg38/sim_{x}')
    if p.exists() and p.is_dir():
        shutil.rmtree(p)



