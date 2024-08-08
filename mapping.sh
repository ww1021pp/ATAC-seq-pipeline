

    mapping
    bowtie2 --mm -p 20 -x $index --very-sensitive -X 2000 -1 ${dir}/$Qc/${sample}_R1.fq.gz -2 $dir/$Qc/${sample}_R2.fq.gz 2> $dir/02_bowtie2_mapping/k0/${sample}.align.log | samtools view -F 4 -Sbhu -@ 20 | samtools sort -@ 20 >$dir/02_bowtie2_mapping/k0/${sample}_mm10.sorted.bam

    2. filter chrM
      samtools idxstats $dir/02_bowtie2_mapping/k0/${sample}_mm10.sorted.bam | cut -f 1 | grep -v  "^chrM" | xargs samtools view $dir/02_bowtie2_mapping/k0/${sample}_mm10.sorted.bam -@ 20 -b> $dir/02_bowtie2_mapping/k0/${sample}_non_mito.bam
    3. filte encode blacklist region
  samtools view $dir/02_bowtie2_mapping/k0/${sample}_non_mito.bam -b -h -o $dir/02_bowtie2_mapping/k0/${sample}_with_blacklsit.bam -U $dir/02_bowtie2_mapping/k0/${sample}_rm_chrM_blck.bam -L /workspace/rsrch2/panpanliu/project/ATAC-seq/Esrrg_ATAC-seq/Scripts/mm10.blacklist.bed
    4. filter low mapping quanlity -q 30
samtools view -h -f 2 -q 30 -@ 20 -u $dir/02_bowtie2_mapping/k0/${sample}_rm_chrM_blck.bam | samtools sort >$dir/02_bowtie2_mapping/k0/${sample}_tmp_filt.bam
    5.  filter duplicate
      picard MarkDuplicates \
      I=$dir/02_bowtie2_mapping/k0/${sample}_tmp_filt.bam \
      O=${dir}/02_bowtie2_mapping/k0/${sample}_marked_duplicates.bam \
      M=${dir}/02_bowtie2_mapping/k0/${sample}_dup.qc \
      VALIDATION_STRINGENCY=LENIENT \
      ASSUME_SORTED=TRUE \
      REMOVE_DUPLICATES=FALSE

     samtools view -F 1024 -@ 20 -b ${dir}/02_bowtie2_mapping/k0/${sample}_marked_duplicates.bam -o ${dir}/02_bowtie2_mapping/k0/${sample}.final.bam
    6. filter bam with insertsize <120
samtools view -h -@ 20 ${dir}/02_bowtie2_mapping/k0/${sample}.final.bam | \
  awk 'substr($0,1,1)=="@" || ($9 <120 && $9 >0) || ($9>=-120 && $9 <0) '| \
  samtools view -b -@ 20 >${dir}/02_bowtie2_mapping/k0/${sample}.filter_is120.bam

7. shift bam for ATAC-seq
alignmentSieve --numberOfProcessors 8 --ATACshift --bam ${dir}/02_bowtie2_mapping/k0/${sample}.filter_is120.bam -o ${dir}/02_bowtie2_mapping/k0/${sample}.tmp_shift.bam 


