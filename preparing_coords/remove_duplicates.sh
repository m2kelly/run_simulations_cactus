#to run on extracted bedfile-to then remove dups before using to call anc etc
#then reintersect with full exon coords
bedfile_in='/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed'
bedfile_out='/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg_nodupcoords.bed'
temp='/home/maria/cactus/to_human_ref/temp'

tail -n +2 $bedfile_in | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' \
> $temp

bedtools subtract -a $temp -b /home/maria/data/segmental_duplicates_cleaned.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3}' | \
bedtools intersect -wa -wb -a - -b $temp \
> $bedfile_out

#have to intersect without flanks, to avoid issues with cut neighbouring regions? 
bedtools subtract -a $temp -b /home/maria/data/segmental_duplicates_cleaned.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3}' | \
bedtools intersect -wa -wb -a - -b $temp | \
awk -F'\t' 'BEGIN{OFS="\t"} ($2>=$5) {print}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3<=$6) {print}' \
> $bedfile_out

rm $temp
'''
#ouput chr start_nodup end_nodup chr start_seq end_seq  hg38 GCA_028858775   Anc4
#problem:sometimes getting start<seq_start
awk -F'\t' 'BEGIN{OFS="\t"} ($6<$3) {print}' '/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg_nodupcoords.bed'
#checking manually
all from adding flanks 
#eg probelm 1: chr1:151,767,220-151,767,230 - from exons with very small gap betwen them 
except for:
chr7-from outputted exon of length 1 in exon_merged_ids - potenially check
chr7    105565446       105565458       chr7    105565275       105565447
otherwise just remove these overlaps
'''

then remove the dup seqs from the bed using removing_duplicates_bed.py
bed_in ='/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg_nodupcoords.bed'
bed_out ='/home/maria/cactus/to_human_ref/coord_filt_nodup_nocpg.bed'
species = ['hg38', 'GCA_028858775','Anc4']


#then reintersect with full exon coords
tail -n +2 /home/maria/cactus/to_human_ref/coord_filt_nodup_nocpg.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3 - 1, $4, $5, $6}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' | \
bedtools intersect -wa -wb -a - -b /home/maria/filter_transcripts/output/exon_merged_ids_sort.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $5, $6, $7, $8 - 1, $9 + 1, $10, $12}' \
> "/home/maria/cactus_target_size/auxillary/extracted_df_nocpg_nodup.bed"



