# Clinical Variant Explorer

Interactive R Shiny application for exploring and filtering annotated genomic variants from multi-sample GRCh38 VCF files. Designed for clinical genetics researchers with minimal coding experience.

## Features

### 📊 Six Interactive Visualization Tabs
1. **Table** — Sortable, filterable variant table with pagination
2. **Manhattan Plot** — Genomic position vs. quality score visualization
3. **Allele Frequency** — MAF distribution histogram + site frequency spectrum
4. **PCA** — Principal component analysis with LD pruning (requires genotype matrix)
5. **Karyotype** — Variant density across chromosomes
6. **IGV** — Integrated genome viewer (click table rows to navigate)

### 🔍 Advanced Filtering
- **Chromosome & region** — Select specific chromosomes or genomic intervals
- **Variant type** — SNP, Indel, MNP
- **Quality scores** — QUAL threshold slider
- **PASS filter** — Show only high-confidence variants
- **Allele frequency** — MAF range selector
- **Gene search** — Autocomplete gene name lookup
- **Consequence** — Filter by SnpEff effect prediction (missense, frameshift, etc.)
- **Impact** — HIGH, MODERATE, LOW, MODIFIER (SnpEff classification)

### 💾 Export Options
- Filtered VCF file
- Variants table (CSV)
- Publication-ready plots (PNG)
- Full PDF report with all visualizations

## Requirements

### Data Assets (must be attached to `/data`)
1. **User VCF + Reference FASTA**
   - Single `.vcf` or `.vcf.gz` file (GRCh38 coordinates)
   - GRCh38 reference FASTA + `.fai` index
   - Optional: `sample_metadata.tsv` with columns `sample_id`, `group`

2. **SnpEff Database** (from `snpeff-grch38-data-builder` capsule)
   - Data Asset ID: `adb864a9-22f9-4c66-bee7-53e573ff10f0`
   - **Mount point**: `/data/snpeff/`
   - Contains: SnpEff v5.2 jar + GRCh38.105 annotation database

### Compute Resources
- **Recommended**: 2x Small (4 cores / 32 GB RAM)
- **Minimum**: Small (2 cores / 16 GB RAM) for smaller VCFs (<10k variants)

## Usage

### 1. Attach Data Assets
Before launching, attach:
- Your VCF + FASTA Data Asset (any mount point)
- `snpeff-grch38-105` Data Asset → **mount to `/data/snpeff/`**

### 2. Launch Cloud Workstation
Click **Start Cloud Workstation** and wait for the Shiny app to launch (~2-5 minutes).

### 3. Access the App
Once running, the app will be accessible at the Cloud Workstation URL on **port 8080**.

### 4. Interact with Visualizations
- Use the left sidebar to apply filters
- Click table rows to navigate IGV to that variant
- Export filtered results via the **Export** button (top-right)

## How It Works

### Stage 1: Preprocessing (`preprocess.R`)
On first launch, the app:
1. Discovers the VCF file in `/data`
2. Runs SnpEff annotation using the mounted GRCh38.105 database
3. Parses the annotated VCF into a tidy data frame
4. Caches the result to `/scratch/variants.parquet` (Parquet format for fast loading)
5. Builds a gene autocomplete index

**Caching**: If you relaunch with the same VCF, preprocessing is skipped (SHA256 check).

### Stage 2: Shiny App (`ui.R`, `server.R`)
- Loads the cached Parquet file
- Applies reactive filters based on sidebar inputs
- Renders six interactive visualization tabs
- Exports filtered data to `/results` on demand

## File Structure
```
/code
  ├── run                    # Launcher script (bash)
  ├── app.R                  # Shiny entrypoint
  ├── preprocess.R           # Stage 1: SnpEff + parsing
  ├── ui.R                   # Shiny UI definition
  ├── server.R               # Shiny server logic
  └── README.md              # This file
/data
  ├── <your_vcf>.vcf.gz      # User-provided VCF
  ├── <reference>.fa         # GRCh38 reference FASTA
  ├── <reference>.fa.fai     # FASTA index
  └── snpeff/                # SnpEff database (mounted)
      ├── snpEff.jar
      ├── snpEff.config
      └── data/GRCh38.105/
/scratch
  ├── variants.parquet       # Cached annotated variants
  ├── source.sha256          # VCF checksum for cache validation
  └── genes.rds              # Gene autocomplete index
/results
  ├── filtered.vcf.gz        # Exported filtered VCF
  ├── filtered_variants.csv  # Exported table
  └── plots/                 # Exported plots
```

## Troubleshooting

### "No VCF file found in /data"
Ensure your VCF Data Asset is attached before launching. Check the Data panel.

### "SnpEff jar not found at /data/snpeff/snpEff.jar"
The SnpEff database Data Asset must be mounted to `/data/snpeff/`. Check the mount path when attaching.

### "Preprocessing is slow"
First-time annotation of large VCFs can take 5-15 minutes. Subsequent launches use the cached Parquet file and start instantly.

### IGV viewer not loading
The current implementation shows a placeholder. Full `igvShiny` integration requires additional configuration (FASTA serving, track loading).

## References
- **SnpEff**: Genomic variant annotation and effect prediction  
  Cingolani et al., *Fly* 2012; DOI: 10.4161/fly.19695

- **karyoploteR**: Chromosome-level visualization in R  
  Gel & Serra, *Bioinformatics* 2017; DOI: 10.1093/bioinformatics/btx346

- **SNPRelate**: PCA and kinship analysis for SNP data  
  Zheng et al., *Bioinformatics* 2012; DOI: 10.1093/bioinformatics/bts606

## License
This capsule is provided as-is for research purposes.
