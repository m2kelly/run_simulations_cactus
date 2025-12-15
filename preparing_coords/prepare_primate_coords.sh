#ran on cluster
#includes cpg islands 
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3 + 1, $4, $6}' '/home/maria/filter_transcripts/output/exon_filt.bed' \
| sort -k1,1V -k2,2n -k3,3n  | uniq \
> /home/maria/cactus/auxilary/coord_flank_filt.bed

#cannot see run using 8 way allignment1!
singularity exec /users/dweghorn/mkelly/containers/cactus_2.9.9_make10.sif \
    cactus-hal2maf ./js_hal2maf_primate \
    /users/dweghorn/mkelly/241_mammals_cactus/data/241-mammalian-2020v2.hal \
    /users/dweghorn/mkelly/cactus_cluster/auxillary/primate.maf.gz \
    --filterGapCausingDupes \
    --binariesMode local \
    --refGenome Homo_sapiens \
    --chunkSize 500000 \
    --batchCores 5 \
    --targetGenomes Homo_sapiens,fullTreeAnc105,fullTreeAnc106,fullTreeAnc107,fullTreeAnc108,fullTreeAnc109,fullTreeAnc110 \
    --batchCount 4 \
    --caching false \
    --filterGapCausingDupes \
    --bedRanges /users/dweghorn/mkelly/cactus_cluster/coord_flank_filt.bed \
    --dupeMode "single" \
    --workDir /no_backup/dweghorn \
    --maxDisk 2000G \
    --batchMemory 4G \
    --restart

#then create bedfile
/home/maria/cactus_target_size/scripts/bed_file_from_my_maf.py
input = '/home/maria/cactus/to_human_ref/primates.maf.gz'
output = "/home/maria/cactus/to_human_ref/primates.bed"
species_output = ['hg38','Anc4','Anc3', 'Anc1', 'Anc0']
ref_species = 'hg38'