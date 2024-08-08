# tumor hep combine call peak
Genrich -t liver12T.sam.Hep_com_NFR,liver34T.sam.Hep_com_NFR -o ../Genrich/tumor_Hep_com.narrowPeak -b ../Genrich/tumor_Hep_com_reads.bed -f ../Genrich/tumor_Hep_com_pq.bed -k ../Genrich/tumor_Hep_com_pileup_p.bed -m 30 -q 0.05 -r -j

# 1. generate common peak
cd /workspace/rsrch1/ychen/Projects/Project01_Human_HCC_sc/cellranger_atac_out/Genrich
# common peak
bedtools intersect -a normal_Hep_com.narrowPeak -b tumor_Hep_com.narrowPeak > ../bedtools/Hep_com_2condition_commonpeak.bed
# unique peak
bedtools subtract -A -a normal_Hep_com.narrowPeak -b tumor_Hep_com.narrowPeak > ../bedtools/normal_Hep_com_unique_peak.bed
bedtools subtract -A -b normal_Hep_com.narrowPeak -a tumor_Hep_com.narrowPeak > ../bedtools/tumor_Hep_com_unique_peak.bed
# merge peak
cat normal_Hep_com.narrowPeak tumor_Hep_com.narrowPeak | sort -k1,1 -k 2,2n > Hep_com_2condition_combined.narrowPeak_sort
bedtools merge -d 10 -i Hep_com_2condition_combined.narrowPeak_sort  > ../bedtools/Hep_com_2condition_merged_peak.bed -k
# Calculate reads count for each merged peaks
## 1. conver bam to bed
cd ~/Project01_Human_HCC_sc/cellranger_atac_out/sam_split_NFR_sort
bedtools bamtobed -i liver12N.sam.Hep_com_NFR > ../bamtobed/liver12N.sam.Hep_com.bed
## 2. calculate reads count for each merged peak
bedtools intersect -a ../bedtools/Hep_com_2condition_merged_peak.bed -b ../bamtobed/
