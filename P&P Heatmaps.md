# Creating sitelist.txt and region_lookup.txt files for each of the probe .bed files
Although there are 2 probe .bed files (probe_ITS2_Hc.bed & probe_ITS_Uni.bed) the script below uses a generic "probe_ITS2.bed"; The workflow is exacly identical for both probe .bed files!

I. probe_ITS2_Hc.bed (Hc probe):
```
hcontortus_chr2_Celeg_TT_arrow_pilon	20855120	20855139
hcontortus_chr2_Celeg_TT_arrow_pilon	20878606	20878625
hcontortus_chr2_Celeg_TT_arrow_pilon	20885517	20885536
hcontortus_chr2_Celeg_TT_arrow_pilon	20892415	20892434
hcontortus_chr2_Celeg_TT_arrow_pilon	20904033	20904052
hcontortus_chr3_Celeg_TT_arrow_pilon	19399240	19399259
hcontortus_chr3_Celeg_TT_arrow_pilon	19406140	19406159
hcontortus_chr3_Celeg_TT_arrow_pilon	19413040	19413059
hcontortus_chr3_Celeg_TT_arrow_pilon	19419940	19419959
hcontortus_chr3_Celeg_TT_arrow_pilon	19426840	19426859
hcontortus_chr3_Celeg_TT_arrow_pilon	19433751	19433770
```

II. probe_ITS_Uni.bed (Univ probe):
```
hcontortus_chr3_Celeg_TT_arrow_pilon	19398931	19398950
hcontortus_chr3_Celeg_TT_arrow_pilon	19405831	19405850
hcontortus_chr3_Celeg_TT_arrow_pilon	19412731	19412750
hcontortus_chr3_Celeg_TT_arrow_pilon	19419631	19419650
hcontortus_chr3_Celeg_TT_arrow_pilon	19426531	19426550
hcontortus_chr3_Celeg_TT_arrow_pilon	19433442	19433461
hcontortus_chr2_Celeg_TT_arrow_pilon	20854811	20854830
hcontortus_chr2_Celeg_TT_arrow_pilon	20878297	20878316
hcontortus_chr2_Celeg_TT_arrow_pilon	20885208	20885227
hcontortus_chr2_Celeg_TT_arrow_pilon	20892106	20892125
hcontortus_chr2_Celeg_TT_arrow_pilon	20903724	20903743
```
```bash
#!/bin/bash -l

#Output folder
output="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/2_convert_bed_table"
#Probe bed file
probe_folder="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data"
#Create it if not existing yet
mkdir -p "$output"

# add a region name column (copy1..copy11; for later sorting)
awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "copy"NR}' "$probe_folder/probe_ITS2.bed" > "$probe_folder/probe_ITS2_named.bed"
# Samtools mpileup will convert from 0 -> +1 coord, so i need to counter that with -1
awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3}' "$probe_folder/probe_ITS2_named.bed" > "$output/sitelist.txt"
# region lookup file generation
awk 'BEGIN{OFS="\t"}{for(p=$2+1;p<=$3;p++) print $1, p, $4}' "$probe_folder/probe_ITS2_named.bed" > "$output/region_lookup.txt"
```
# Creating mpileups
```bash
#!/bin/bash -l

set -euo pipefail

module load samtools/1.20

# def dirs
bam_files="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/6_mark_duplicate"
site_file="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/2_convert_bed_table"
output="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/3_pileup_probe_bases"
refer_genome='/cfs/klemming/projects/supr/naiss2025-23-66/variants_explore_2025-2/raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

mkdir -p "$output"

for i in "$bam_files"/*.bam; do
    base_name=$(basename "$i" .bam)

    samtools mpileup -f "$refer_genome" -l "$site_file/sitelist.txt" -aa "$i" \
        > "$output/${base_name}.pileup.txt"

    echo "pileup done for $base_name"
done
```
# Parsing pileup.txt(s) and assembling nucleotide counts in a table (A, C, T, G, N and del)
```python
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
lookup = pd.read_csv("region_lookup.txt", sep="\t", names=["chrom","pos","copy"])

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
# Combining each .tsv into 1 df & calculating fraction % of non reference alleles
```r
library(tidyverse)

setwd("C:/Users/pauli/Desktop/SLU/data")

# Generate a list metafile
files<-list.files(pattern = "\\.raw_counts\\.tsv")

combined <- files %>%
  set_names() %>% ## every vector is a vector with it's own name
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "source_file") %>% ## looping through every file applying a function to bind rows of each into a big table (_dfr data frame) ## adds a column .id = "source file" to id where that file came from
  mutate(population = str_remove(source_file, "\\.sort\\.markdup\\.raw_counts.tsv$")) # edit the column source file
```
```r
# relative position for each copy (nts 1 to 20)
## assign unique groups of chrom and copy so that operations (min calculation) are done on pos (taking into account chrom and copy combo) instead of just on pos alone (not accounting for chrom or copy).

combined <- combined %>%
  group_by(chrom, copy) %>%
  mutate(relpos=pos-min(pos) +1) %>%
  ungroup()
```
```r
# CHECK: Is any ref base across the copies different or are all identical? More than 1 distinct ref allele per copies?

ref_check <- combined %>%
  group_by(relpos) %>%
  summarise(n_distinct_ref = n_distinct(ref)) %>%
  filter(n_distinct_ref > 1)
```
```r
# Combine mpileup reads for every relative position across copies for all pops (version 2 = N and deletions)
agg<- combined %>%
  group_by(population, relpos) %>%
  summarise(
    ref = first(ref),
    count_A =  sum(count_A),
    count_C = sum(count_C),
    count_G = sum(count_G),
    count_T = sum(count_T),
    count_N = sum(count_N),
    count_del = sum(count_del),
    .groups= "drop"
  ) %>%
  mutate(called_bases = count_A + count_C + count_G + count_T + count_N + count_del)
```
```r
# Compute fractions

agg <- agg %>%
  rowwise() %>% # perform below operations rowwise
  mutate(
    count_ref = get(paste0("count_", toupper(ref))), # count_ref variable, which is count_X, depending on what ref value is in that row
    frac_nonref = (called_bases - count_ref) / called_bases # creates another column where called bases - counted ref alleles / called bases = frac of non-ref alleles
  ) %>%
  ungroup()
```
```r
# Remove "ampure" tags
agg$population<- gsub("-ampure.*","", agg$population)
```
```r
# Plotting the data table

heatmap<-ggplot(agg, aes(x=relpos, y= population, fill= frac_nonref)) +
  geom_tile(color="white", linewidth = 0.2)+
  
  scale_fill_gradient2(low="#3B4CC0", mid="#F7F7F7", high = "#B40426", midpoint = 0.5, limits = c(0,1),
                       name="Fraction\nnon-reference alleles") +
  scale_x_continuous(breaks =1:20)+
  labs(x = "Relative *HC* probe 5' position within a copy of rDNA gene (1-20 nt)",
       y="Population",
       title = NULL) +

theme_minimal() +
  theme(
    axis.text.x= element_text(size = 8),
    axis.text.y= element_text(size= 8),
    axis.title.x = element_markdown()
  )
```
# Creating sitelist.txt and region_lookup.txt files for each of the primer .bed files

Once again, there are 2 primer .bed files (primers_ITS2_Hc.bed & primers_ITS_Uni.bed), the script below uses a generic "primer_sites_named.bed"; The workflow is exacly identical for both primer .bed files!

I. primers_ITS2_Hc.bed (Hc primers):
```
hcontortus_chr3_Celeg_TT_arrow_pilon	19399145	19399167	copy1	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19406045	19406067	copy2	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19412945	19412967	copy3	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19419845	19419867	copy4	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19426745	19426767	copy5	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19433656	19433678	copy6	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20855025	20855047	copy7	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20878511	20878533	copy8	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20885422	20885444	copy9	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20892321	20892343	copy10	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20903939	20903961	copy11	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19399292	19399311	copy1	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19406192	19406211	copy2	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19413092	19413111	copy3	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19419992	19420011	copy4	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19426892	19426911	copy5	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19433803	19433822	copy6	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20855172	20855191	copy7	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20878658	20878677	copy8	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20885569	20885588	copy9	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20892467	20892486	copy10	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20904085	20904104	copy11	Reverse
```

II. primers_ITS_Uni.bed (Univ primers):
```
hcontortus_chr3_Celeg_TT_arrow_pilon	19398881	19398901	copy1	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19405781	19405801	copy2	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19412681	19412701	copy3	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19419581	19419601	copy4	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19426481	19426501	copy5	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19433392	19433412	copy6	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20854761	20854781	copy7	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20878247	20878267	copy8	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20885158	20885178	copy9	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20892056	20892076	copy10	Forward
hcontortus_chr2_Celeg_TT_arrow_pilon	20903674	20903694	copy11	Forward
hcontortus_chr3_Celeg_TT_arrow_pilon	19398972	19398989	copy1	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19405872	19405889	copy2	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19412772	19412789	copy3	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19419672	19419689	copy4	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19426572	19426589	copy5	Reverse
hcontortus_chr3_Celeg_TT_arrow_pilon	19433483	19433500	copy6	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20854852	20854869	copy7	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20878338	20878355	copy8	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20885249	20885266	copy9	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20892147	20892164	copy10	Reverse
hcontortus_chr2_Celeg_TT_arrow_pilon	20903765	20903782	copy11	Reverse
```
```bash
awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3}' "primer_sites_named.bed" > "primer_sitelist.txt"
```

```bash
awk 'BEGIN{OFS="\t"}{for(p=$2;p<=$3;p++) print $1, p, $4, $5}' "primer_sites_named.bed" \
    > "primer_region_lookup.txt"
```
# Creating mpileups containing
```bash
#!/bin/bash -l

set -euo pipefail

module load samtools/1.20

# def dirs
bam_files="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/6_mark_duplicate"
site_file="/cfs/klemming/home/p/pabs0001/naiss2025-23-66/paulius_ITS2_primers/data/2.5_convert_bed_table_primersHC"
output="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/4_pileup_primer_bases"
refer_genome='/cfs/klemming/projects/supr/naiss2025-23-66/variants_explore_2025-2/raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

mkdir -p "$output"

for i in "$bam_files"/*.bam; do
    base_name=$(basename "$i" .bam)

    samtools mpileup -f "$refer_genome" -l "$site_file/primer_sitelist.txt" -aa "$i" \
        > "$output/${base_name}.pileup.txt"

    echo "pileup done for $base_name"
done

```
# Parsing pileup.txt(s) and assembling nucleotide counts in a table (A, C, T, G, N and del)
```python
# Set dir wih files
import os
os.chdir(r"C:\Users\pauli\Desktop\SLU\data\primers_HC") 

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
lookup = pd.read_csv("primer_region_lookup.txt", sep="\t", names=["chrom","pos","copy","direction"])

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
    merged = merged[["chrom","pos","copy","direction","ref","depth",
                      "count_A","count_C","count_G","count_T","count_N","count_del"]]
    # save it!
    merged.to_csv(f"{pop_name}.primer_raw_counts.tsv", sep="\t", index=False)
    print(f"done: {pop_name}")
```

# Combining each .tsv into 1 df & calculating fraction % of non reference alleles

```r
library(tidyverse)

setwd("C:/Users/pauli/Desktop/SLU/data/primers_HC/")

# Generate a list metafile
files<-list.files(pattern = "\\_raw_counts\\.tsv")

combined2 <- files %>%
  set_names() %>% ## every vector is a vector with it's own name
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "source_file") %>% ## looping through every file applying a function to bind rows of each into a big table (_dfr data frame) ## adds a column .id = "source file" to id where that file came from
  mutate(population = str_remove(source_file, "\\.sort\\.markdup\\.primer_raw_counts.tsv$")) # edit the column source file
```
```r
 Forward and Reverse should both start at 5', rather than at the ''left/upstream most" base
combined2 <- combined2 %>%
  group_by(chrom, copy, direction) %>%
  mutate(
    relpos = if_else(direction == "Reverse",
                      max(pos) - pos + 1,   # flip numbering for reverse-strand primers
                      pos - min(pos) + 1)   # normal left-to-right numbering for forward
  ) %>%
  ungroup()
```
```r
# CHECK: Is any ref base across the copies different or are all identical (grouped by direction as well)? More than 1 distinct ref allele per copies?

ref_check <- combined2 %>%
  group_by(direction, relpos) %>%
  summarise(n_distinct_ref = n_distinct(ref), .groups = "drop") %>%
  filter(n_distinct_ref > 1)
```
```r
# Combine mpileup reads for every relative position across copies for all pops (version 2 = N and deletions)
agg2<- combined2 %>%
  group_by(population, direction, relpos) %>%
  summarise(
    ref = first(ref),
    count_A =  sum(count_A),
    count_C = sum(count_C),
    count_G = sum(count_G),
    count_T = sum(count_T),
    count_N = sum(count_N),
    count_del = sum(count_del),
    .groups= "drop"
  ) %>%
  mutate(called_bases = count_A + count_C + count_G + count_T + count_N + count_del)
```
```
# Compute fractions

agg2 <- agg2 %>%
  rowwise() %>% # perform below operations rowwise
  mutate(
    count_ref = get(paste0("count_", toupper(ref))), # count_ref variable, which is count_X, depending on what ref value is in that row
    frac_nonref = (called_bases - count_ref) / called_bases # creates another column where called bases - counted ref alleles / called bases = frac of non-ref alleles
  ) %>%
  ungroup()
```
```r
# Remove "ampure" tags

agg2$population<- gsub("-ampure.*","", agg2$population)
```
```r
library(patchwork)

heatmap2<-ggplot(agg2, aes(x=relpos, y= population, fill= frac_nonref)) +
  geom_tile(color="white", linewidth = 0.2)+
  
  scale_fill_gradient2(low="#3B4CC0", mid="#F7F7F7", high = "#B40426", midpoint = 0.5, limits = c(0,1),
                       name="Fraction\nnon-reference alleles") +
  scale_x_continuous(breaks =1:max(agg2$relpos))+
  labs(x = "Relative *HC* primer 5' position within a copy of rDNA gene (1-20 nt R; 1-23 nt F)",
       y="Population",
       title = NULL) +

facet_wrap(~direction, ncol=1)+
  
theme_minimal() +
  theme(
    axis.text.x= element_text(size = 8),
    axis.text.y= element_text(size= 8),
     axis.title.x = element_markdown()
  )

heatmap2
```
# Combining heatmap plots
```
#Combine primer and probe plots into 1 thru patchwork
comb_plot<-heatmap+heatmap2 +
  #plot_annotation(title = 'Non-reference allele frequencies in the HC-specific primer and probe binding sites across the different populations') +
  plot_layout(guides = "collect", axes="collect")
comb_plot
```
