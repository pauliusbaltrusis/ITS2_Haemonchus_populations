# Add taxids
```bash
sed 's/^>/>kraken:taxid|6289|/'  6289_Haemonchus_contortus.fna         > tagged_6289.fna
sed 's/^>/>kraken:taxid|45464|/' 45464_Teladorsagia_circumcincta.fna   > tagged_45464.fna
sed 's/^>/>kraken:taxid|6319|/'  6319_Trichostrongylus_colubriformis.fna > tagged_6319.fna
sed 's/^>/>kraken:taxid|40349|/' 40349_Trichostrongylus_axei.fna      > tagged_40349.fna
sed 's/^>/>kraken:taxid|27828|/' 27828_Cooperia_oncophora.fna         > tagged_27828.fna
sed 's/^>/>kraken:taxid|63234|/' 63234_Oesophagostomum_venulosum.fna  > tagged_63234.fna
sed 's/^>/>kraken:taxid|6317|/'  6317_Ostertagia_ostertagi.fna        > tagged_6317.fna
sed 's/^>/>kraken:taxid|9940|/'  9940_Ovis_aries.fna > tagged_9940.fna
```
# Build Kraken2 DB
```bash
set -euo pipefail
shopt -s nullglob

module load PDC/26.03
module load kraken2/2.17.1

genome_folder="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/nem_genomes"
working_folder="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/7_Kraken2_output"

cd $working_folder

DBNAME=GIN_nemabiome_db

kraken2-build --download-taxonomy --skip-maps --db "$DBNAME" ### can run this from the login node with --use-ftp flag

for f in "$genome_folder"/tagged_*.fna; do
  kraken2-build --add-to-library "$f" --db "$DBNAME"
done

kraken2-build --build --db "$DBNAME" --threads 8
```

# Run Kraken2
```bash
set -euo pipefail
shopt -s nullglob

module load PDC/26.03
module load kraken2/2.17.1

reads_folder="/cfs/klemming/home/p/pabs0001/naiss2025-23-66/paulius_ITS2_primers/2_trim_fastp"
working_folder="/cfs/klemming/projects/supr/naiss2025-23-66/paulius_ITS2_primers/data/7_Kraken2_output"

cd $working_folder

DBNAME=GIN_nemabiome_db

for r1 in "$reads_folder"/*_R1.fastq.gz; do
  r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
  sample=$(basename "$r1" _fastp_trim_R1.fastq.gz)
  
  kraken2 --db "$DBNAME" \
    --threads 8 \
    --paired \
    --gzip-compressed \
    --report "${sample}.kreport" \
    --output "${sample}.kraken2" \
    --use-names \
    --confidence 0.0 \
    "$r1" "$r2"
done
```
# Bracken
```bash
set -euo pipefail

module load PDC/26.03
module load bioinfo-tools
module load kraken2/2.17.1
module load Bracken/2.6.2

reports_folder="/cfs/klemming/home/p/pabs0001/naiss2025-23-66/paulius_ITS2_primers/data/7_Kraken2_output"

cd $reports_folder

DBNAME=GIN_nemabiome_db
READLEN=150

bracken-build -d "$DBNAME" -t 8 -k 35 -l 150

## per sample using generated .kreport files
for kreport in *.kreport; do
  sample=$(basename "$kreport" .kreport)
  bracken -d "$DBNAME" -i "$kreport" -o "${sample}.bracken" -l S -r $READLEN
done
```
# Plotting
```r
library(tidyverse)


bracken_dir<-"C:/Users/pauli/Desktop/SLU/data/Kraken2 db/0_confidence reports"

# Generate a list metafile
files_bracken<-list.files(bracken_dir, pattern = "\\.bracken$", full.names = TRUE)

just_sample_names<-basename(files_bracken) %>% str_remove("\\.bracken$")



bracken_all <- map_dfr(set_names(files_bracken, just_sample_names),
                       ~read_tsv(.x, show_col_types = FALSE), .id= "sample")


bracken_tidy <- bracken_all %>%
  rename(species=name) %>%
  mutate(species = str_trim(species)) %>%
  select(sample,species, taxonomy_id, new_est_reads, fraction_total_reads)

# Normalize to 100% before adding in the unclassified reads

bracken_prct_norm <- bracken_tidy %>%
  group_by(sample) %>%
  mutate(prct= new_est_reads / sum(new_est_reads) * 100) %>%
  ungroup()

# Grabbing the unclassified reads from kraken2 report "U"

kreport_dir <- "C:/Users/pauli/Desktop/SLU/data/Kraken2 db/0_confidence reports/kraken2 reports"

# extract U reads from each report using a function
get_unclassified_prct <- function(sample) {
  kreport_file <- file.path(kreport_dir, paste0(sample, ".kreport"))
  kreport<- read_tsv(kreport_file, col_names = c( "prct", "reads_clade", "reads_direct", "rank", "taxid", "name"), show_col_types = FALSE)
  u_row <-kreport %>% filter(rank=="U")
  if (nrow(u_row)==0) return(0)
  u_row$prct[1]
}

unclassified_prct<-tibble(sample=just_sample_names) %>%
  mutate(unclass_prct=map_dbl(sample, get_unclassified_prct))

# Removing the unclassfied read % from the total per sample

bracken_prct_with_unclassified <- bracken_prct_norm %>%
  left_join(unclassified_prct, by ="sample") %>%
  mutate(prct = prct * (100-unclass_prct)/100) %>%
  select(sample, species, taxonomy_id, new_est_reads, prct)


unclassified_row <- unclassified_prct %>%
  transmute(sample, species = "Unclassified", taxonomy_id =NA, new_est_reads=NA, prct=unclass_prct)

# Joining the U% row with bracken_pcrt_with_unclassified so every adds up to 100%

bracken_prct_with_unclassified <- bind_rows(bracken_prct_with_unclassified, unclassified_row)

# Remove "ampure" tags

bracken_prct_with_unclassified$sample<- gsub("-ampure.*","", bracken_prct_with_unclassified$sample)

#Plotting
species_order<- plot_data %>%
  filter(species!="Unclassified") %>%
  distinct(species) %>%
  arrange(species) %>%
  pull(species) %>%
  as.character()

n_species<-length(species_order)


species_palette <- c("#3B4252", "#5E81AC", "#EBCB8B", "#8FBCBB",
                       "#A3BE8C", "#B48EAD", "#D08770", "#BF616A")[1:n_species]

stopifnot(length(species_palette) ==n_species) # will error if wrong

plot_data <- plot_data %>%
  mutate(species=factor(species, levels= c(species_order, "Unclassified")))


species_colors<- setNames(c(species_palette, "grey50"),
                          c(species_order, "Unclassified"))

stopifnot(!anyNA(species_colors)) # will error if wrong

p_stacked<-ggplot(plot_data, aes(x=sample, y=prct, fill=species)) +
  geom_col(width=0.7) +
  scale_fill_manual(values= species_colors)+
  scale_y_continuous(labels=scales::percent_format(scale=1), expand=c(0,0))+
  labs(x=NULL,
       y="Relative abundance(%)",
       fill="Species",
       title="") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.x=element_blank())
p_stacked

p_facet<- ggplot(plot_data, aes(x=sample, y=prct, fill=species)) +
  geom_col(show.legend = FALSE)+
  facet_wrap(~species, scales="free_y")+
  scale_fill_manual(values= species_colors)+
  labs(
    x=NULL,
    y="Relative abundance (%)",
    title=""
  ) +
  theme_minimal(base_size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        strip.text = element_text(face="bold"))
  
p_facet
```
