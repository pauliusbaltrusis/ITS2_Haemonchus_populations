# Generating a .bed file (=amplicon_regions.bed) 

(containing the loci for Hc ITS2 amplicon (167bp) chr, start, end and copy #)
```bash
hcontortus_chr3_Celeg_TT_arrow_pilon	19399145	19399311	copy1
hcontortus_chr3_Celeg_TT_arrow_pilon	19406045	19406211	copy2
hcontortus_chr3_Celeg_TT_arrow_pilon	19412945	19413111	copy3
hcontortus_chr3_Celeg_TT_arrow_pilon	19419845	19420011	copy4
hcontortus_chr3_Celeg_TT_arrow_pilon	19426745	19426911	copy5
hcontortus_chr3_Celeg_TT_arrow_pilon	19433656	19433822	copy6
hcontortus_chr2_Celeg_TT_arrow_pilon	20855025	20855191	copy7
hcontortus_chr2_Celeg_TT_arrow_pilon	20878511	20878677	copy8
hcontortus_chr2_Celeg_TT_arrow_pilon	20885422	20885588	copy9
hcontortus_chr2_Celeg_TT_arrow_pilon	20892321	20892486	copy10
hcontortus_chr2_Celeg_TT_arrow_pilon	20903939	20904104	copy11
```
# Sitelist and region_lookup files

```bash
awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3}' "$site_file/amplicon_regions.bed" > "$site_file/amplicon_sitelist.txt"
awk 'BEGIN{OFS="\t"}{for(p=$2;p<=$3;p++) print $1, p, $4}' "$site_file/amplicon_regions.bed" > "$site_file/amplicon_region_lookup.txt"
```
# Generating mpileups

```bash
module load samtools/1.20

for bam in "$bam_files"/*.bam; do
    base_name=$(basename "$bam" .bam)
    samtools mpileup -f "$refer_genome" -l "$site_file/amplicon_sitelist.txt" -aa "$bam" \
        > "$output/${base_name}.amplicon.pileup.txt"
    echo "done: $base_name"
done
```
# Parsing mpileups into A, C, T, G, N, del counts and aggregating across copies
```python
# Set dir wih files
import os
os.chdir(r"C:\Users\pauli\Desktop\SLU\data\HC_amplicon")

#### Version 2 A C T G N INDELS 

# Import packages
import glob
import pandas as pd
import re

def parse_pileup_base_string(bases, ref):
    bases = re.sub(r'\^.', '', bases)
    bases = bases.replace('$', '')
    bases = re.sub(r'[+-](\d+)([ACGTNacgtn]+)', '', bases)
    counts = {"A":0, "C":0, "G":0, "T":0, "N":0, "del":0}
    for b in bases:
        if b in '.,':
            if ref.upper() in counts:
                counts[ref.upper()] += 1
        elif b in '*#': # Deletions have their own labelling but not insertions. Insertions (such as +2AG, exist outside the reference labelling entirely =>
            # and it's not possible to assign then to particular positions i.e. 1-20 nts for example)
            counts["del"] += 1
        elif b.upper() == 'N':
            counts["N"] += 1
        elif b.upper() in counts:
            counts[b.upper()] += 1
    return counts

# Use a look up file containing all bases in the range for 11 ITS2 copies
lookup = pd.read_csv("amplicon_region_lookup.txt", sep="\t", names=["chrom","pos","copy"])

# Loop through .pileup(s) converting every line into a string of matches/mismatches to ref and recording them in a table
for pileup_file in glob.glob("*.pileup.txt"):
    pop_name = pileup_file.replace(".pileup.txt", "")
    rows = []
    with open(pileup_file) as f:
        for line in f:
            # remove trailing new lines and replace with tabs
            fields = line.rstrip("\n").split("\t")
            # all 4 columns are pulled out and pos and depth are converted to integers
            chrom, pos, ref, depth = fields[0], int(fields[1]), fields[2], int(fields[3])
            # bases column is 4th and need to check it even exists (not empty with depth 0), otherwise its empty
            bases = fields[4] if len(fields) > 4 and depth != "0" else ""
            # calling the previous function to count bases
            counts = parse_pileup_base_string(bases, ref)
            # adding the date to the df
            rows.append([chrom, pos, ref, depth, counts["A"], counts["C"], counts["G"], counts["T"],
                         counts["N"], counts["del"]])

    df = pd.DataFrame(rows, columns=["chrom","pos","ref","depth",
                                       "count_A","count_C","count_G","count_T",
                                       "count_N","count_del"])
    # Merging lookup file with every base for 11 copies with the generated df 
    merged = df.merge(lookup, on=["chrom","pos"], how="left")
    merged = merged[["chrom","pos","copy","ref","depth",
                      "count_A","count_C","count_G","count_T","count_N","count_del"]]
    # save it!
    merged.to_csv(f"{pop_name}.raw_counts.tsv", sep="\t", index=False)
    print(f"done: {pop_name}")
```
```r
library(tidyverse)
library(ggtext)

setwd("C:/Users/pauli/Desktop/SLU/data/HC_amplicon")

# Generate a list metafile
files_ampl<-list.files(pattern = "\\.raw_counts\\.tsv")

combined_ampl <- files_ampl %>%
  set_names() %>% ## every vector is a vector with it's own name
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "source_file") %>% ## looping through every file applying a function to bind rows of each into a big table (_dfr data frame) ## adds a column .id = "source file" to id where that file came from
  mutate(population = str_remove(source_file, "\\.sort\\.markdup\\.raw_counts.tsv$")) # edit the column source file

```
```r
# relative position for each copy (nts 1 to 1xx)
## assign unique groups of chrom and copy so that operations (min calculation) are done on pos (taking into account chrom and copy combo) instead of just on pos alone (not accounting for chrom or copy).

combined_ampl <- combined_ampl %>%
  group_by(chrom, copy) %>%
  mutate(relpos=pos-min(pos) +1) %>%
  ungroup()
  
  # push the copy10 and 11 frame +1 past relpos 45 to account for the deletion in 45th position in the amplicon
  combined_ampl <- combined_ampl %>%
  mutate(
    relpos = case_when(
      copy %in% c("copy10", "copy11") & relpos > 45 ~ relpos + 1,
      TRUE ~ relpos
    )
  )

  # 45th base is T as expected, so deletion is between 45-47 nt
combined_ampl %>% filter(relpos==45) %>% select(ref)
  
```
```r
# Exclude the four positions where there are multiple alleles
positions_to_exclude <- c(40, 62, 68, 145)

# remove the rows containing positions to exclude
combined_ampl_without <- combined_ampl %>%
  filter(!relpos %in% positions_to_exclude) ### some positions excluded due to low-level inter-copy sequence variation
```
```r
agg_ampl<- combined_ampl_without %>%
  group_by(population, relpos) %>%
  summarise(
    ref = dplyr::first(ref),
    count_A =  sum(count_A),
    count_C = sum(count_C),
    count_G = sum(count_G),
    count_T = sum(count_T),
    count_N = sum(count_N),
    count_del = sum(count_del),
    .groups= "drop"
  ) %>%
  mutate(called_bases = count_A + count_C + count_G + count_T + count_N + count_del,
         freq_A = count_A / called_bases,
         freq_C = count_C / called_bases,
         freq_G = count_G / called_bases,
         freq_T = count_T / called_bases,
         freq_N = count_N / called_bases,
         freq_del = count_del / called_bases)
```
Comments: I Deletion 45-47 meant i needed to push the frame by +1 after relpos >45. II Pos (40, 62, 68, 145) were excluded due to naturally-occurring inter-copy variation at those reference sites
#Building population x feature matrix for PCA and running PCA
```r
library(tidyr)

wide <- agg %>%
  select(population, relpos, freq_A, freq_C, freq_G, freq_T) %>%
  pivot_longer(cols = starts_with("freq_"), names_to = "base", values_to = "freq") %>%
  mutate(feature = paste0("pos", relpos, "_", base)) %>%
  select(population, feature, freq) %>%
  pivot_wider(names_from = feature, values_from = freq)

pca_matrix <- as.matrix(wide[, -1])
rownames(pca_matrix) <- wide$population
```
```r
# remove any columns with zero variance (e.g. perfectly conserved positions across all populations)
pca_matrix_filtered <- pca_matrix[, apply(pca_matrix, 2, var) > 0]

pca_result <- prcomp(pca_matrix_filtered, scale. = TRUE)

summary(pca_result)   # variance explained by each PC
```
# Plotting
## PCA
```r
library(ggplot2)
library(RColorBrewer)
# convert to a df and return pop names in rownames
pca_df<- as.data.frame(pca_result$x)
pca_df$population<- rownames(pca_df)


pca_plot<-ggplot(pca_df, aes(x=PC1, y=PC2, color=population))+
  geom_point(size=4)+
  scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(length(unique(pca_df$population))))+
  labs(x=paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100,1), "%)"),
  y=paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100,1), "%)"),
  title="") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_line(color="grey40")
  )

pca_plot
  ```
## Depth at every Hc amplicon position per population
```r
# Supplementary - plot the depth at every pos (1-167 for all populations) - barplot histogram

depth_agg<- agg_ampl %>% select(population, relpos, ref, called_bases) %>%
  mutate(population=gsub("-ampure.*","",population)) %>% mutate(population=gsub(".sort.markdup.*","",population))
#rename rows old <- new
depth_agg$population<-real_names[match(depth_agg$population, temp_names)]

# order
depth_agg$population <- factor(depth_agg$population, levels = ranked_names)

# plot
depth_plot<- ggplot(depth_agg, aes(x=relpos, y=called_bases, color=population)) +
  geom_col(show.legend = FALSE)+
  facet_wrap(~population, scales="free_y")+
  scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(length(unique(depth_agg$population))))+
  labs(
    x=NULL,
    y="Depth",
    title=""
  ) +
  theme_minimal(base_size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        strip.text = element_text(face="bold"))
  
depth_plot
```
