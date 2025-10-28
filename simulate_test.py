from NeutralEvolutionSimulator import NeutralEvolutionSimulator

for x in range(5,10):
    print(f'running simulation {x}')
    simulator = NeutralEvolutionSimulator(
        n_subs=72047,
        start_speci_file='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/processed_gene_seqs.pkl',
        cbase_output='/home/maria/run_simulations_cactus/auxiliary/cbase.csv',
        syn_probs_dir='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/syn_target',
        non_syn_probs_dir='/home/maria/run_simulations_cactus/auxiliary/anc4_HCLCA/nonsyn_target',
        output=f'/home/maria/run_simulations_cactus/output_anc4_hg38/simulated_genes_{x}.pkl',
        gene_name_map = '/home/maria/filter_transcripts/output/gene_name_id.csv')
    simulator.load_cbase_data()