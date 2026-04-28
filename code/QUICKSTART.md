# Quick Start Guide — Clinical Variant Explorer

## Prerequisites

Before launching this capsule, you need **two Data Assets** attached:

### 1. Your VCF + Reference Data
Create a Data Asset containing:
- Your VCF file (`.vcf` or `.vcf.gz`) with GRCh38 coordinates
- GRCh38 reference FASTA (`.fa`) + index (`.fai`)
- (Optional) `sample_metadata.tsv` with columns: `sample_id`, `group`

**Attach to**: Any mount point (e.g., `/data/vcf_data`)

### 2. SnpEff Database (already attached to this capsule)
- Data Asset ID: `adb864a9-22f9-4c66-bee7-53e573ff10f0`
- Name: `snpeff-grch38-105`
- **Mount point**: `/data/snpeff/` ✅ (already configured)

---

## Step-by-Step Launch

### Step 1: Attach Your VCF Data Asset
1. Go to the **Data** panel in the capsule editor
2. Click **+ Attach Data Asset**
3. Select your VCF Data Asset
4. Mount to any location (e.g., `vcf_input`)
5. Click **Save**

### Step 2: Validate Setup (optional but recommended)
Open a Cloud Workstation terminal and run:
```bash
Rscript /code/validate_setup.R
```

This checks:
- ✓ VCF file exists
- ✓ SnpEff database is accessible
- ✓ Java is installed

### Step 3: Launch the App
**Option A: Cloud Workstation (Recommended)**
1. Click **Start Cloud Workstation**
2. Wait 2-5 minutes for environment setup
3. The Shiny app will launch automatically on port 8080
4. Access via the Cloud Workstation URL

**Option B: Reproducible Run (batch mode)**
- Use this if you only want to preprocess and export filtered data
- Results will be written to `/results`
- No interactive UI; requires command-line arguments

### Step 4: Interact with the App
Once the app loads:
- **Sidebar**: Apply filters (chromosome, gene, quality, etc.)
- **Table tab**: Browse variants; click rows to navigate IGV
- **Other tabs**: Visualizations update automatically based on filters
- **Export button**: Download filtered VCF, CSV, or plots

---

## Expected Behavior

### First Launch (~5-15 minutes)
- **Stage 1**: SnpEff annotates your VCF (~3-10 min for large files)
- **Stage 2**: VCF is parsed and cached to `/scratch/variants.parquet`
- **Stage 3**: Shiny app launches and loads cached data

### Subsequent Launches (~30 seconds)
- Cache is reused (same VCF detected via SHA256)
- Shiny app starts immediately

---

## Troubleshooting

### "No VCF file found in /data"
**Solution**: Attach your VCF Data Asset before launching.

### "SnpEff jar not found"
**Solution**: The SnpEff database must be mounted to `/data/snpeff/`. Check the Data panel and verify mount path.

### App takes too long to load
- First-time annotation is slow (10-15 min for WGS-scale VCFs)
- Check Cloud Workstation logs for progress
- Consider downsampling your VCF for testing

### Out of memory
**Solution**: Increase compute to **Medium** (8 cores / 64 GB RAM) for very large VCFs (>100k variants).

---

## Output Files (written to `/results`)

When you click **Export** in the app:
- `filtered.vcf.gz` — VCF containing only variants that pass your filters
- `filtered_variants.csv` — Table export
- `plots/` — PNG exports of each visualization
- `variant_report.pdf` — Full report with all plots (if selected)

---

## Technical Notes

### Reproducibility
To ensure reproducibility:
1. **Release this capsule** after configuring it with your data
2. Attach the same Data Assets every time
3. Run via **Reproducible Run** (not Cloud Workstation) for batch exports

### Performance
- Small VCF (<10k variants): X-Small compute (1 core / 8 GB)
- Medium VCF (10k-100k variants): Small compute (2 cores / 16 GB)
- Large VCF (>100k variants): 2x Small or Medium (4+ cores / 32+ GB)

### Cache Invalidation
To force re-annotation (e.g., if you update the VCF):
1. Delete `/scratch/variants.parquet` and `/scratch/source.sha256`
2. Relaunch the app

---

## Need Help?

Check the main README.md for detailed documentation:
```bash
cat /code/README.md
```

Or run the validation script:
```bash
Rscript /code/validate_setup.R
```

---

## Example Workflow

1. Upload your cohort VCF + reference FASTA to Code Ocean as a Data Asset
2. Attach the Data Asset to this capsule
3. Launch Cloud Workstation
4. Wait for preprocessing to complete
5. Filter variants by gene (e.g., *BRCA1*), impact (HIGH), and MAF (<0.01)
6. Export filtered VCF for further analysis
7. Share the filtered results or publish the capsule with a Release

---

**Ready to explore your variants!** 🧬
