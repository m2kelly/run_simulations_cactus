import pandas as pd

'''
seq_file ='/home/maria/run_simulations_cactus/auxillary_anc3/anc3_nodup_nocpg_alligned.bed'
output_file = '/home/maria/run_simulations_cactus/auxillary_anc3/anc3_nodup_nocpg_flank.bed'
'''

seq_file ='/home/maria/run_simulations_cactus/auxillary_wg/anc4_nodup_nocpg_alligned.bed'
output_file = '/home/maria/run_simulations_cactus/auxillary_wg/anc4_nodup_nocpg.bed'

#ensure output empty
with open(output_file,'w') as f:
    f.write('')

#header = ['chr', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy','start', 'end', 'gene','gene_name','strand' ]
header = ['chr','start', 'end', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy','exon_start', 'exon_end', 'gene','gene_name','strand']
chunk_iter = pd.read_csv(seq_file, chunksize=10000, sep='\t', names=header)

for chunk in chunk_iter:
    chunk_dfs =[]
    #add flanks
    #chunk.drop(columns=['chr_copy','gene','gene_name'],inplace=True)
    #print problems 
    
    chunk['start'] = chunk['start'].astype(int)
    chunk['seq_start'] = chunk['seq_start'].astype(int)
    chunk['end'] = chunk['end'].astype(int)
    #add flanks
    chunk['start'] -=1
    chunk['end'] +=1
    chunk['exon_start'] = chunk['exon_start'].astype(int)
    chunk['exon_end'] = chunk['exon_end'].astype(int)
    chunk['exon_start'] -= 1
    chunk['exon_end'] += 1

    no_left_mask = chunk['seq_start']>chunk['start']
    no_right_mask  = chunk['seq_end']<chunk['end']
    
    
    chunk_no_flank = chunk[no_left_mask|no_right_mask].copy()
    if not chunk_no_flank.empty:
        chunk_no_flank['seq'] = chunk_no_flank.apply(lambda row: 'F'+row['seq'][row['start']+1-row['seq_start']: row['end']-1-row['seq_start']]+'F',axis=1)
        chunk_no_flank['target_seq']=chunk_no_flank.apply(lambda row:'F'+ row['target_seq'][row['start']+1-row['seq_start']: row['end']-1-row['seq_start']]+'F',axis=1)
        chunk_dfs.append(chunk_no_flank)
        chunk = chunk[~(no_left_mask|no_right_mask)].copy()
        
        

        no_left_mask = chunk['seq_start']>chunk['start']
        no_right_mask  = chunk['seq_end']<chunk['end']

   
    #then add F's instead of flanks 
    chunk_no_flank_right = chunk[no_right_mask].copy()
    if not chunk_no_flank_right.empty:
        chunk_no_flank_right['seq'] = chunk_no_flank_right.apply(lambda row: row['seq'][row['start']-row['seq_start']: row['end']-1-row['seq_start']]+'F',axis=1)
        chunk_no_flank_right['target_seq']=chunk_no_flank_right.apply(lambda row: row['target_seq'][row['start']-row['seq_start']: row['end']-1-row['seq_start']]+'F',axis=1)
        chunk_dfs.append(chunk_no_flank_right)
        chunk = chunk[~no_right_mask].copy()
        no_left_mask = chunk['seq_start']>chunk['start']
        no_right_mask  = chunk['seq_end']<chunk['end']
    

    chunk_no_flank_left = chunk[no_left_mask].copy()
    if not chunk_no_flank_left.empty:
        chunk_no_flank_left['seq'] = chunk_no_flank_left.apply(lambda row: 'F' + row['seq'][row['start'] + 1 -row['seq_start']: row['end']-row['seq_start']],axis=1)
        chunk_no_flank_left['target_seq']=chunk_no_flank_left.apply(lambda row: 'F' + row['target_seq'][row['start'] + 1 -row['seq_start']: row['end']-row['seq_start']],axis=1)
        chunk_dfs.append(chunk_no_flank_left)
        chunk = chunk[~no_left_mask].copy()

    #extract seqs
    if not chunk.empty:
        print(f'percent with flanks = {len(chunk)/100}')
        chunk['seq'] = chunk.apply(lambda row: row['seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
        chunk['target_seq']=chunk.apply(lambda row: row['target_seq'][row['start']-row['seq_start']: row['end']-row['seq_start']],axis=1)
        chunk_dfs.append(chunk)
    full_chunk = pd.concat(chunk_dfs,axis=0)
        
    #add flanks to exon starts 
    
    full_chunk[[ 'chr','start', 'end', 'target_seq','seq', 'chr_copy','exon_start', 'exon_end', 'gene','strand']].to_csv(output_file,header=None, mode='a',sep='\t',index=False)

