#anc3
'''
input=/home/maria/effective_pop/data/anc4_anc3.bed
output='/home/maria/run_simulations_cactus/auxillary_anc3/anc3_nodup_nocpg_alligned.bed'
coords='/home/maria/run_simulations_cactus/auxillary_anc3/exons_nodup_nocpg_flank_coords.bed'
exon_coords='/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
exon_nooverlap='/home/maria/filter_transcripts/output/exon_filt.bed'
'''

input=/home/maria/effective_pop/data/HCLCA_wg.bed
output='/home/maria/run_simulations_cactus/auxillary_wg/anc4_nodup_nocpg_alligned.bed'
coords='/home/maria/run_simulations_cactus/auxillary_wg/exons_nodup_nocpg_flank_coords.bed'
exon_coords='/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
exon_nooverlap='/home/maria/filter_transcripts/output/exon_filt.bed'

tail -n +2 $input | \
cut -f1-3 | \
bedtools intersect -a - -b $exon_nooverlap | \
cut -f1-3 | \
bedtools subtract -a - -b /home/maria/data/segmental_duplicates_cleaned.bed | \
bedtools subtract -a - -b /home/maria/data/cpgIslandExt.bed \
> $coords

tail -n +2 $input | \
bedtools intersect -wa -wb -a $coords -b - | \
cut -f 1-3,5-8 | \
bedtools intersect -wa -wb -a - -b $exon_coords | \
uniq \
> $output

#target=anc4
#seq=anc3
#[ 'chr','extract start', 'extract end', 'seq_start','seq_end', 'target_seq','seq', 'chr_copy','exon start', 'exon end', 'gene','gene_name','strand']


#test


input=/home/maria/effective_pop/data/HCLCA_wg.bed
output='/home/maria/run_simulations_cactus/auxillary_wg/anc4_nodup_nocpg_alligned.bed'
coords='/home/maria/run_simulations_cactus/auxillary_wg/exons_nodup_nocpg_flank_coords.bed'
exon_coords='/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
exon_nooverlap='/home/maria/filter_transcripts/output/exon_filt.bed'

tail -n +2 $input | \
cut -f1-3 | \
bedtools intersect -a - -b $exon_nooverlap | \
cut -f1-3 | \
bedtools subtract -a - -b /home/maria/data/segmental_duplicates_cleaned.bed | \
bedtools subtract -a - -b /home/maria/data/cpgIslandExt.bed | \
awk '{len = $3 - $2; if (len > 0) s += len} END{print s+0}' 

awk '{len = $3 - $2; if (len > 0) s += len} END{print s+0}' $exon_coords