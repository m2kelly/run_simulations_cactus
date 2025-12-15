from concurrent.futures import ProcessPoolExecutor
from GeneContentCounter import GeneContentCounter

def run_for_sim(sim):
    input_file = f'/home/maria/run_simulations_cactus/output_anc4_hg38/simulated_genes_{sim}.pkl'
    output_file = f'/home/maria/run_simulations_cactus/output_anc4_hg38/sim_{sim}_heuristics'
    generator = GeneContentCounter(input_file,output_file)
    generator.run()


with ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(run_for_sim, list(range(1,100)))