sample <- c('liver12N', 'liver12T',  'liver34N', 'liver34T')

data <- c()
for (x in sample){
  data1 <- read.table(paste0('readscount/', x, '_2condition_merged_peak.count'))
  names(data1) <- c('chr', 'start', 'end', 'count')
  data = data.frame(
    id = paste(data1$chr, data1$start, data1$end, sep = '.'),
    chr = data1$chr,
    start =  data1$start,
    end =  data1$end,
    count = data1$count,
    sample = x) %>%
    rbind(., data)
}

data_convert <- spread(data, sample, count) %>% column_to_rownames(var = 'id')
head(data_convert)
        
colData <- data.frame(sample = sample, condition = c('Normal', 'Tumor', 'Normal', 'Tumor')) %>% column_to_rownames(var = 'sample')
colData

# perform DESeq2 --------
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data_convert[,4:7], # expression matrix
                              colData = colData, # metadata
                              design = ~ condition) # compare between condition
# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
# run DESeq2
dds <- DESeq(dds)
# Check the coefficients for the comparison
resultsNames(dds)
# Generate results object
res <- as.data.frame(results(dds, name = "condition_Tumor_vs_Normal"))
res <- merge(data_convert, res, by = 'row.names')
res$summit_up200 <- res$start + ceiling((res$end - res$start)/2)-250
res$summit_down200 <-  res$start + ceiling((res$end - res$start)/2)+250
ggplot(data = res, aes(x = log2FoldChange, y = -log2(padj))) + geom_point() + theme_bw() + geom_hline(yintercept = -log(0.01))
write.table(res, './diff_peak/hep_com_Tumor_vs_Normal_peak.txt', sep = '\t', quote = F, row.names = F)

head(res)
res$rpkm_liver12N <- res$liver12N * 10^6/(sum(res$liver12N) * (res$end - res$start + 1))
res$rpkm_liver12T <- res$liver12T * 10^6/(sum(res$liver12T) * (res$end - res$start + 1))
res$rpkm_liver34N <- res$liver34N * 10^6/(sum(res$liver34N) * (res$end - res$start + 1))
res$rpkm_liver34T <- res$liver34T * 10^6/(sum(res$liver34T) * (res$end - res$start + 1))
write.table(res, './diff_peak/hep_com_Tumor_vs_Normal_peak_with_rpkm.txt', sep = '\t', quote = F, row.names = F)

tumor_hep_com_up <- read.table('/workspace/rsrch1/ychen/Projects/Project01_Human_HCC_sc/cellranger_atac_out/bedtools/tumor_Hep_com_unique_peak.bed')
tumor_hep_com_up <- tumor_hep_com_up[,1:3]
names(tumor_hep_com_up) <- c('chr', 'start', 'end')
tumor_hep_com_up <- tumor_hep_com_up %>%
  mutate(summit_up200 = start + ceiling((end - start)/2)-250) %>%
  mutate(summit_down200 = start + ceiling((end - start)/2)+250) %>%
  mutate(id = paste(chr, start, end, sep = '.')) %>%
  dplyr::select(chr, summit_up200, summit_down200, id)
write.table(tumor_hep_com_up, './diff_peak/hep_com_Tumor_unique_peak_501bp.txt',sep = '\t', quote = F, row.names = F, col.names = F)

tumor_hep_com_up <- read.table('/workspace/rsrch1/ychen/Projects/Project01_Human_HCC_sc/cellranger_atac_out/bedtools/normal_Hep_com_unique_peak.bed')
tumor_hep_com_up <- tumor_hep_com_up[,1:3]
names(tumor_hep_com_up) <- c('chr', 'start', 'end')
tumor_hep_com_up <- tumor_hep_com_up %>%
  mutate(summit_up200 = start + ceiling((end - start)/2)-250) %>%
  mutate(summit_down200 = start + ceiling((end - start)/2)+250) %>%
  mutate(id = paste(chr, start, end, sep = '.')) %>%
  dplyr::select(chr, summit_up200, summit_down200, id)
write.table(tumor_hep_com_up, './diff_peak/hep_com_Normal_unique_peak_501bp.txt',sep = '\t', quote = F, row.names = F, col.names = F)
