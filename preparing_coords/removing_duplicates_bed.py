import pandas as pd

bed_in ='/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg_nodupcoords.bed'
bed_out ='/home/maria/cactus/to_human_ref/coord_filt_nodup_nocpg.bed'
species = ['hg38', 'GCA_028858775','Anc4']

def cut_seq(start,end, seq_start,seq_end, seq_list):
    cut_start = start-seq_start
    cut_end = end-seq_start
    if end > seq_end:
        raise ValueError('end > seq end')
    if start < seq_start:
        raise ValueError('start < seq start')
    seqs_cut =[]
    for seq in seq_list:
        seqs_cut.append(seq[cut_start:cut_end])
    return seqs_cut

df = pd.read_csv(bed_in, names=['chr','start','end','chr_copy','seq_start','seq_end']+species,sep='\t')
mask = (df['seq_start']!=df['start']) | (df['seq_end']!=df['end'])
df_to_cut = df[mask].copy()

df_to_cut[species] = df_to_cut.apply(lambda row: pd.Series(cut_seq(row['start'],row['end'], row['seq_start'],row['seq_end'], row[species])),axis=1)
df_keep = df[~mask].copy()

df_filtered = pd.concat([df_keep,df_to_cut])
df_filtered[['chr','start','end']+species].to_csv(bed_out, index=False, sep='\t')