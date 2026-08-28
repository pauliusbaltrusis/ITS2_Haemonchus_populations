# Add taxids
sed 's/^>/>kraken:taxid|6289|/'  6289_Haemonchus_contortus.fna         > tagged_6289.fna
sed 's/^>/>kraken:taxid|45464|/' 45464_Teladorsagia_circumcincta.fna   > tagged_45464.fna
sed 's/^>/>kraken:taxid|6319|/'  6319_Trichostrongylus_colubriformis.fna > tagged_6319.fna
sed 's/^>/>kraken:taxid|40349|/' 40349_Trichostrongylus_axei.fna      > tagged_40349.fna
sed 's/^>/>kraken:taxid|27828|/' 27828_Cooperia_oncophora.fna         > tagged_27828.fna
sed 's/^>/>kraken:taxid|63234|/' 63234_Oesophagostomum_venulosum.fna  > tagged_63234.fna
sed 's/^>/>kraken:taxid|6317|/'  6317_Ostertagia_ostertagi.fna        > tagged_6317.fna
sed 's/^>/>kraken:taxid|9940|/'  9940_Ovis_aries.fna > tagged_9940.fna

# Build Kraken2 DB

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

# Run Kraken2

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

# Bracken

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

# Plotting
```
