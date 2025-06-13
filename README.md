# STRUCTURE Analysis Data Preparation Pipeline

This pipeline processes filtered SNP genotype data and prepares it for analysis with the software **[STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html)**. The workflow includes marker extraction, SNP filtering, genotype cleaning, allele conversion, data transposition, sample filtering, and final formatting into STRUCTURE-compatible format.

---

## Directory Structure

```
structure_data_prep_pipeline/
├── Snakefile
├── config.yaml
├── data/
│   ├── txt_file_with_labs_5PercAF_5PercH_20PercmissInSNPandSample_StepSetAs8_E2ref.txt  # Input file containing SNPs with metadata
│   ├── unwanted_samples.txt         # List of sample IDs to exclude
├── scripts/
│   ├── filter_snps.py               # Filters biallelic SNPs
│   ├── convert_to_real_alleles.py   # Converts 'A'/'B' to REF/ALT alleles
│   ├── format_structure.py          # Final numeric encoding for STRUCTURE
├── results/
│   └── structure_174300_snps.txt    # Final STRUCTURE-ready file
```

---

## Dependencies

- **Snakemake** ≥ 7.0
- **Python** ≥ 3.6
- **pandas**, **datamash**, **grep**, **sed**, **cut**, **paste**

---

## How to Run

```bash
snakemake --cores 4
```

---

## Analysis Steps

### 1. **Raw Data Cleanup**

**Rule:** `clean_raw_data`  
**Tool:** `sed`, `cut`  
**Description:** Removes quotes, extracts relevant SNP marker fields, and saves to `data_with_markers.txt`.

---

### 2. **Filter Valid SNPs**

**Rule:** `filter_snps`  
**Script:** `scripts/filter_snps.py`  
**Description:** Filters for SNPs where both REF and ALT alleles are single nucleotides. Keeps the header.

---

### 3. **Replace Heterozygous and Missing Codes**

**Rule:** `replace_missing_codes`  
**Tool:** `sed`, `cut`, `paste`  
**Description:**

- Replaces ambiguous codes:
  - H → `-`, D → `-`, C → `-`
- Splits metadata and genotype fields
- Merges into `snp_matrix_cleaned.txt`

---

### 4. **Convert A/B Genotypes to Real Alleles**

**Rule:** `convert_to_real_alleles`  
**Script:** `scripts/convert_to_real_alleles.py`  
**Description:** Replaces `A` with reference allele and `B` with alternate allele for each SNP.

---

### 5. **Transpose and Filter Unwanted Samples**

**Rule:** `transpose_and_filter`  
**Tool:** `cut`, `datamash`, `grep`  
**Description:**

- Removes the REF/ALT columns, keeping only genotype matrix.
- Transposes rows → columns.
- Removes samples listed in `unwanted_samples.txt`.

---

### 6. **Numeric Encoding for STRUCTURE**

**Rule:** `format_structure`  
**Script:** `scripts/format_structure.py`  
**Description:** Converts alleles to STRUCTURE numeric codes:

- Three-allele sites: [missing = `-9`, major = `1`, minor = `0`]
- Two-allele sites: [major = `1`, minor = `0`]

---

## Output

- `results/structure_174300_snps.txt`: Final STRUCTURE-ready SNP matrix with samples as rows and SNPs as columns.

---

## Notes

- This pipeline assumes SNPs are encoded as `A`, `B`, or other characters.
- Missing/ambiguous genotypes (`H`, `D`, `C`) are masked using `-`.
- Adapt the unwanted sample list in `data/unwanted_samples.txt` as needed.

---

# AMOVA and FST Analyses

This project contains R scripts and supporting files for performing genetic diversity analyses, including **heterozygosity statistics**, **FST calculations**, and **AMOVA** (Analysis of Molecular Variance), using SNP genotype data of various populations.

## Directory Structure

```
amova_fst_analyses/
├── data/
│ ├── pop_file.txt # Population metadata (individual IDs and population labels)
│ └── 174300_r_snps.txt # SNP genotype matrix used in Structure analysis
│
├── results/ # Output directory for all analysis results (created by script)
│ ├── genind_summary.txt
│ ├── Ho_values.csv
│ ├── He_values.csv
│ ├── overall_stats.csv
│ ├── He_Ho_distribution.pdf
│ ├── fst_wc.txt
│ ├── pairwise_fst_wc84.csv
│ ├── pairwise_fst_nei87.csv
│ ├── amova_result.txt
│ ├── amova_randtest.txt
│ └── amova_randtest_plot.pdf
│
└── scripts/
└── genetic_diversity_analysis.R # Main R script performing FST and AMOVA analyses
```

## Analyses Performed

Replace missing data in a structure-format file:
`sed 's/\b-9\b/NA/g' structure_data_prep_pipeline/results/structure_174300_snps.txt > amova_fst_analyses/data/174300_r_snps.txt`

The script `genetic_diversity_analysis.R` executes the following steps:

1. **Setup & Data Loading**

   - Sets the working directory to the script location.
   - Loads genotype data and population metadata from the `../data/` folder.

2. **Data Preparation**

   - Formats SNP genotype matrix for analysis.
   - Converts genotype data into `genind` and `hierfstat` objects using the `adegenet` and `hierfstat` packages.

3. **Genetic Diversity Statistics**

   - Computes basic diversity metrics:
     - **Observed heterozygosity (Ho)**
     - **Expected heterozygosity (He)**
     - **Overall statistics**
   - Plots histograms of Ho and He distributions.

4. **FST Estimations**

   - Calculates:
     - Global FST using the **Weir and Cockerham** method.
     - **Pairwise FST** values using both **WC84** and **Nei87** methods.

5. **AMOVA Analysis**

   - Converts genotype data to `genclone` format and applies **clone correction**.
   - Performs **AMOVA** using population groupings.
   - Conducts a **permutation test** (n = 1000) to assess significance of AMOVA results.
   - Plots the null distribution and test statistic from the permutation test.

6. **Output**
   - All results (statistics, tables, and plots) are saved in the `../results/` directory.

## Requirements

The R script uses the following R packages:

- `hierfstat`
- `adegenet`
- `dplyr`
- `poppr`
- `ggplot2`

Ensure all packages are installed before running the script:

```r
install.packages(c("hierfstat", "adegenet", "dplyr", "poppr", "ggplot2"))
```
