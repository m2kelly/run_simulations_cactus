#get coords without cpgs islands and with flanks:
bedtools subtract -a '/home/maria/filter_transcripts/output/exon_filt.bed' -b /home/maria/data/cpgIslandExt.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $6}'  \
| sort -k1,1V -k2,2n -k3,3n  | uniq \
> /home/maria/cactus/auxilary/coord_flank_filt_nocpg.bed

#extract sequences using /home/maria/cactus_container/prepare.sh:
cactus-hal2maf ./js_hal2maf8_filt_nocpg /mnt/8_way_primate/8-t2t-apes-2023v2.hal /mnt/to_human_ref/coord_filt_nocpg.maf.gz \
    --filterGapCausingDupes \
    --refGenome hg38 \
    --chunkSize 100000 \
    --batchCores 20 \
    --targetGenomes Anc4,GCA_028858775.2\
    --batchCount 4 \
    --caching false \
    --filterGapCausingDupes \
    --bedRanges /mnt/auxilary/coord_flank_filt_nocpg.bed \
    --dupeMode "single"

#extract bed file from maf
run /home/maria/cactus_target_size/scripts/bed_file_from_my_maf.py 
with:
input = '/home/maria/cactus/to_human_ref/coord_filt_nocpg.maf.gz'
output = "/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed2"
species_output = ['hg38','Anc4','GCA_028858775']
ref_species = 'hg38'

#no such coord_filt_no_ins_nocpg.bed2, assume renamed to coord_filt_no_ins_nocpg.bed2


#then reintersect with full exon coords
tail -n +2 /home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3 - 1, $4, $5, $6}' | \
awk -F'\t' 'BEGIN{OFS="\t"} ($3>$2) {print}' | \
bedtools intersect -wa -wb -a - -b /home/maria/filter_transcripts/output/exon_merged_ids_sort.bed | \
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $5, $6, $7, $8 - 1, $9 + 1, $10, $12}' \
> "/home/maria/cactus_target_size/auxillary/extracted_df2_nocpg.bed"
