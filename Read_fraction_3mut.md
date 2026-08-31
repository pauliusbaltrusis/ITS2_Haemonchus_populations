# Estimating read fraction % containing 3 or more mutated alleles (compared to the ref genome) in the loci overlapping the Hc & Univ probes
## Padding the region by 150 bp on each side, extracting smaller .bam files, indexed .bai and ref genome

probe_regions.bed is either "probe_ITS2_Hc.bed" or "probe_ITS_Uni.bed" from [P&P Heatmaps.md](https://github.com/pauliusbaltrusis/ITS2_Haemonchus_populations/blob/main/P%26P%20Heatmaps.md)

```bash
module load samtools/1.20 

# pad by ~150bp so reads starting/ending just outside the window still get pulled in 

awk 'BEGIN{OFS="\t"}{print $1, $2-150, $3+150}' "$site_file/probe_regions.bed" > "$site_file/probe_regions_padded.bed" 

extract_dir="/cfs/klemming/.../data/6_probe_snp_extract" mkdir -p "$extract_dir" 

for bam in "$bam_files"/*.bam; do 
	base_name=$(basename "$bam" .bam) samtools view -b -h -L "$site_file/probe_regions_padded.bed" "$bam" \ 
	> "$extract_dir/${base_name}.probe_region.bam" 
samtools index "$extract_dir/${base_name}.probe_region.bam" 
echo "extracted: $base_name" done 
#check .bam file sizes
du -sh "$extract_dir"/*.bam # check sizes before downloading

```
## Function to count the mismatches for every read (hits; and retain if min_mismatches>=3) overlapping probe loci

```r
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)

# Load reference genome (needs a .fai index alongside it - run `samtools faidx setwd('C:/Users/pauli/Desktop/SLU/data/probe_HC/padded_regions_filter')

ref_fasta <- FaFile('C:/Users/pauli/Desktop/SLU/data/HC reg genome indexed/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa')

regions<-read_tsv('C:/Users/pauli/Desktop/SLU/data/probe_HC/probe_ITS2.bed', col_names = c('chrom', 'start', 'end')) %>% mutate(copy=paste0("copy", 1:11))

count_mismatches_per_read <- function(bam_path, chrom, start, end, min_mismatches=3) {
  
  region_gr <- GRanges(chrom, IRanges(start, end))
  
  # pull reads overlapping the region where the probe binds, including seq + cigar
  
  param<- ScanBamParam(which=region_gr, what=c("qname", "seq"))
  alns <- readGAlignments(bam_path, param= param, use.names=TRUE)
  
  if (length(alns) == 0) {
    return(data.frame(total_reads=0, reads_with_hit=0))
    
  }
    #ref seq for the exact window
    ref_seq<-as.character(getSeq(ref_fasta, region_gr))
    
    # seq is stored as a column in mcols(alns)
    seqs <- mcols(alns)$seq    
    
    #project every reads seq onto ref coords
    # insertions get dropped ( no ref position)
    # Deletions become "-" chars (so they align 1-1 with ref)
    layered <- sequenceLayer(seqs, cigar(alns), from="query", to= "reference", D.letter="-", N.letter = "-")
    
    # stack all reads into a matrix spanning start : end
    # padding with "+" when a read doesnt cover that position
    stacked<- stackStrings(layered, from=start, to=end,
                           shift = start(alns)-1, Lpadding.letter="+", Rpadding.letter = "+")
    stacked_mat<-as.matrix(stacked) # rows = reads, cols = position in a region
    ref_vec <- strsplit(ref_seq, "")[[1]]
    
    total_reads <-nrow(stacked_mat)
    
    hits<-0
    
    for (i in seq_len(total_reads)) {
      read_row <- stacked_mat[i, ]
      covered<- read_row!= "+" #ignore positions this read does not span
      mismatches <- sum(toupper(read_row[covered]) != toupper(ref_vec[covered]))
      # "-" deletiojn will already count as mismatch, because it's not a ref base
      if (mismatches >=min_mismatches) hits <- hits +1
    }
    
    data.frame(total_reads=total_reads, reads_with_hit = hits)  
    }
    
```
## Loop the function through all smaller/truncated .bam files

```r
bam_files <- list.files(pattern = "\\.probe_region\\.bam$")

results <- list()
for (bam_path in bam_files) {
  pop_name <- sub("\\.probe_region\\.bam$", "", bam_path)

  for (i in seq_len(nrow(regions))) {
    res <- count_mismatches_per_read(
      bam_path, regions$chrom[i], regions$start[i], regions$end[i], min_mismatches = 3
    )
    results[[length(results) + 1]] <- data.frame(
      population = pop_name,
      copy = regions$copy[i],
      chrom = regions$chrom[i],
      start = regions$start[i],
      end = regions$end[i],
      total_reads = res$total_reads,
      reads_with_ge3_mismatches = res$reads_with_hit,
      pct = ifelse(res$total_reads > 0, 100 * res$reads_with_hit / res$total_reads, NA)
    )
  }
  cat("done:", pop_name, "\n")
}

results_df <- do.call(rbind, results)
write.table(results_df, "probe_region_ge3snp_reads.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```
## Collapsing over all 11 copies
```r
library(dplyr)

pop_summary <- results_df %>% group_by(population) %>%
 summarise( total_reads = sum(total_reads), reads_with_ge3_mismatches = sum(reads_with_ge3_mismatches) ) %>%
 mutate(pct = 100 * reads_with_ge3_mismatches / total_reads)

write.table(pop_summary, "probe_region_ge3snp_reads_per_population.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```
## Plotting
```r
 Merge pop_summary tables
pop_summary_Hc<- pop_summary %>% mutate(pct_hc=pct) %>% select(population,pct_hc)
pop_summary_Univ<-pop_summary_uni %>% mutate(pct_univ=pct) %>% select(population, pct_univ)

pop_summaries<-inner_join(pop_summary_Hc, pop_summary_Univ, by="population")

pop_summaries$population<-  gsub("\\-ampure.*", "", pop_summaries$population)
pop_summaries$population<-  gsub("\\.sort\\.markdup", "", pop_summaries$population)

temp_names <- c("Sample_XJ-4126-1", "Sample_XJ-4126-2", "Sample_XJ-4126-3",
                   "Sample_XJ-4126-4", "Sample_XJ-4126-5", "Sample_XJ-4126-6",
                   "Sample_XJ-4126-7", "Sample_XJ-4126-8", "Sample_XJ-4126-9",
                   "Sample_XJ-4126-10", "Sample_XJ-4126-11", "Sample_XJ-4126-12",
                   "Sample_XJ-4126-13", "Sample_XJ-4126-14", "Sample_XJ-4126-15",
                   "Sample_XJ-4126-16", "Sample_XJ-4126-17")  #  temp names


real_names<-c("A2", "F1", "C", "D", "E", "A1", "G1", "B", "G2", "F2", "H", "N", "I", "L", "K", "M", "J")


ranked_names<-c("A1","A2","B",'C','D','E','F1','F2','G1','G2','H','I','J','K','L',"M","N")


pop_summaries$population<-real_names[match(pop_summaries$population, temp_names)]

pop_summaries$population <- factor(pop_summaries$population, levels = ranked_names)

pop_summaries_longer<-pivot_longer(pop_summaries, names_to = "group", values_to = "prct", cols = starts_with("pct"))

# Plotting
pop_summary_plot<-ggplot(pop_summaries_longer, aes(x=population, y=prct, fill= group))+
  geom_col(position=position_dodge(width=0.8), width = 0.7, alpha=0.7) +
  labs(x= "", y="Reads with >=3 mutations %", fill="")+
  scale_fill_manual(labels = c("Hc", "Univ"), values = c("#3B4CC0", "#06402B"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1))

pop_summary_plot
```

