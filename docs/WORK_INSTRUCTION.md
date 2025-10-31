# Krocus Work Instruction

## Purpose

This document provides step-by-step instructions for using Krocus to perform multi-locus sequence typing (MLST) from long-read sequencing data.

## Scope

This work instruction applies to:
- PacBio sequencing data
- Oxford Nanopore sequencing data
- Any uncorrected long-read FASTQ files

## Prerequisites

### Software Requirements
- Python 3.3 or higher
- Krocus installed (`pip3 install krocus`)
- Internet access (for database download)

### Input Requirements
- FASTQ file with long reads (can be gzipped)
- Minimum 10x genome coverage recommended
- Read length: 100bp minimum

### Knowledge Requirements
- Basic command-line usage
- Understanding of MLST concepts
- Familiarity with your bacterial species

## Safety and Quality

### Data Quality Checks

Before analysis, verify:
1. FASTQ file is not corrupted (can be opened/unzipped)
2. File size is reasonable (>1MB for typical bacterial genome)
3. Read quality scores present (not all N's)

### Important Notes

⚠️ **Warnings:**
- Very low quality reads (<Q7 average) may not produce results
- Mixed samples will give unpredictable results
- Incorrect species database will give incorrect results

## Procedure

### Step 1: Database Setup

#### 1.1 List Available Species

```bash
krocus_database_downloader --list_species
```

**Expected Output:** List of bacterial species with MLST schemes

**Troubleshooting:**
- If no output: Check internet connection
- If error: Verify krocus is installed correctly

#### 1.2 Download Species Database

```bash
krocus_database_downloader \
    --species "Salmonella enterica" \
    --output_directory Salmonella_db \
    --verbose
```

**Parameters:**
- `--species`: Exact species name from list (case-sensitive)
- `--output_directory`: Directory to store database
- `--verbose`: Show download progress

**Expected Output:**
```
Downloading "https://pubmlst.org/..." ... done
Alleles found:
['aroC', 'dnaN', 'hemD', 'hisD', 'purE', 'sucA', 'thrA']
```

**Success Criteria:**
- Directory created with `.tfa` files
- `profile.txt` file present
- No error messages

**Troubleshooting:**
- "Species not found": Check spelling, use exact name from list
- Download errors: Check internet, try again later
- Permission denied: Check directory write permissions

### Step 2: Quality Check

#### 2.1 Verify Database

```bash
ls -lh Salmonella_db/
```

**Expected Output:**
- Multiple `.tfa` files (one per gene)
- One `profile.txt` file
- File sizes: 10KB - 500KB typical

#### 2.2 Verify FASTQ File

```bash
# Check file exists and size
ls -lh reads.fastq.gz

# View first few lines (if uncompressed)
head -n 8 reads.fastq

# Or for gzipped
gunzip -c reads.fastq.gz | head -n 8
```

**Expected Output:**
```
@read_id_1
ATCGATCGATCG...
+
IIIIIIIIII...
```

**Quality Check:**
- File size >1MB for typical analysis
- Reads present (not empty file)
- Quality scores present (fourth line each read)

### Step 3: Run Krocus

#### 3.1 Basic Analysis

```bash
krocus Salmonella_db/ reads.fastq.gz
```

**Expected Runtime:**
- Small dataset (<100MB): 1-5 minutes
- Medium dataset (100-500MB): 5-15 minutes
- Large dataset (>500MB): 15-60 minutes

**Expected Output:**
```
10      45.23   aroC(1)*  dnaN(5)  hemD(12)  hisD(1)  purE(2)*  sucA(1)  thrA(3)
323     97.23   aroC(10)  dnaN(5)  hemD(12)  hisD(1)  purE(2)   sucA(1)  thrA(3)
```

**Output Format:**
- Column 1: Sequence type (ST)
- Column 2: K-mer coverage percentage
- Column 3+: Gene(allele) pairs

**Symbols:**
- `*` = Partial match (some k-mers missing)
- No symbol = Complete match

#### 3.2 Advanced Options

For higher quality or different read types, adjust k-mer size:

```bash
# High quality reads (error rate <5%)
krocus Salmonella_db/ reads.fastq.gz --kmer 15

# Medium quality reads (error rate ~10%) - DEFAULT
krocus Salmonella_db/ reads.fastq.gz --kmer 11

# Low quality reads (error rate >15%)
krocus Salmonella_db/ reads.fastq.gz --kmer 7
```

### Step 4: Save Results

#### 4.1 Save to File

```bash
krocus Salmonella_db/ reads.fastq.gz \
    --output_file results.txt
```

**Verification:**
```bash
cat results.txt
```

#### 4.2 Save Filtered Reads (Optional)

Save reads that match MLST genes for further analysis:

```bash
krocus Salmonella_db/ reads.fastq.gz \
    --output_file results.txt \
    --filtered_reads_file matching_reads.fastq
```

**Use Cases:**
- De novo assembly of specific genes
- Quality assessment
- Further analysis

### Step 5: Interpret Results

#### 5.1 Good Result Example

```
323     97.23   aroC(10)  dnaN(5)  hemD(12)  hisD(1)  purE(2)  sucA(1)  thrA(3)
```

**Interpretation:**
- ✅ High coverage (97.23%)
- ✅ Clear ST identified (323)
- ✅ No partial matches
- **Conclusion**: High confidence result

#### 5.2 Partial Match Example

```
10      78.45   aroC(1)*  dnaN(5)  hemD(12)*  hisD(1)  purE(2)*  sucA(1)  thrA(3)
```

**Interpretation:**
- ⚠️ Lower coverage (78.45%)
- ⚠️ Multiple partial matches (*)
- ⚠️ ST may be correct but uncertain

**Actions:**
1. Check k-mer size (may be too large)
2. Check read quality
3. Check coverage depth
4. If >50% partial matches, consider increasing read depth

#### 5.3 No Result or Very Low Coverage

```
ND      12.34   aroC(1)*  dnaN(?)*  hemD(?)*  ...
```

**Interpretation:**
- ❌ Very low coverage
- ❌ ST not determined (ND)
- ❌ Result not reliable

**Troubleshooting:**
1. Reduce k-mer size: `--kmer 7`
2. Check if correct species database used
3. Verify adequate sequencing coverage
4. Check for contamination

### Step 6: Real-time Analysis (Optional)

For real-time monitoring during sequencing:

```bash
tail -f active_run.fastq | krocus Salmonella_db/ - --print_interval 10
```

**Options:**
- `-` reads from stdin
- `--print_interval 10` updates every 10 reads

**Stop Condition:**
When coverage plateaus or ST stabilizes, analysis is complete.

### Step 7: Documentation

Record the following in lab notebook:

1. **Sample Information:**
   - Sample ID
   - Date of analysis
   - Sequencing platform used
   - Database version/date downloaded

2. **Analysis Parameters:**
   - K-mer size used
   - Krocus version (`krocus --version`)
   - Command line used

3. **Results:**
   - Sequence type identified
   - Coverage percentage
   - Any partial matches
   - Interpretation/confidence level

4. **Files Generated:**
   - Output file location
   - Filtered reads file (if applicable)
   - Database directory used

## Quality Control

### Acceptance Criteria

Results are acceptable if:
- ✅ Coverage >80%
- ✅ ST clearly identified (not ND)
- ✅ <3 genes with partial matches

Results are questionable if:
- ⚠️ Coverage 50-80%
- ⚠️ Multiple partial matches
- ⚠️ ST changes significantly during real-time analysis

Results should be rejected if:
- ❌ Coverage <50%
- ❌ ST = ND (not determined)
- ❌ All or most genes show partial matches

### Validation

For critical samples, validate results by:
1. Running with different k-mer sizes
2. Comparing with reference method (if available)
3. Checking against known samples
4. Confirming with Sanger sequencing of MLST genes

## Troubleshooting Guide

### Problem: No Output

**Causes:**
- Incorrect database path
- No matching reads
- Empty FASTQ file

**Solutions:**
1. Verify database directory exists and contains files
2. Check FASTQ file is not empty
3. Try reducing `--min_fasta_hits` to 5
4. Try reducing k-mer size

### Problem: Very Low Coverage

**Causes:**
- Wrong species database
- K-mer size too large for read quality
- Insufficient sequencing depth
- Sample contamination

**Solutions:**
1. Verify correct species database
2. Reduce k-mer size: `--kmer 7`
3. Check sequencing depth (use FastQC or similar)
4. Check for mixed samples

### Problem: Many Partial Matches

**Causes:**
- K-mer size too large
- Novel alleles
- Read errors in specific regions

**Solutions:**
1. Reduce k-mer size
2. Accept result if coverage >90% overall
3. Consider Sanger sequencing for confirmation

### Problem: Analysis Too Slow

**Causes:**
- Very large file
- Low-performance computer

**Solutions:**
1. Use sampling: Extract subset of reads first
2. Increase `--print_interval` to reduce output overhead
3. Use streaming from compressed file:
   ```bash
   gunzip -c reads.fastq.gz | krocus db/ -
   ```

## Appendix A: Command Reference

### Quick Commands

```bash
# List species
krocus_database_downloader -l

# Download database
krocus_database_downloader -s "Species name" -o output_dir

# Basic analysis
krocus database_dir/ reads.fastq

# Save results
krocus database_dir/ reads.fastq -o results.txt

# Real-time analysis
tail -f reads.fastq | krocus database_dir/ - -p 10

# Adjust for low quality
krocus database_dir/ reads.fastq --kmer 7

# Adjust for high quality
krocus database_dir/ reads.fastq --kmer 15
```

## Appendix B: File Formats

### Database Directory Structure

```
database_dir/
├── profile.txt          # ST definitions
├── gene1.tfa            # Allele sequences
├── gene2.tfa
└── ...
```

### profile.txt Format

```
ST    gene1    gene2    gene3    ...
1     1        1        1        ...
2     1        1        2        ...
```

### Output Format

```
ST    Coverage    gene1(allele)  gene2(allele)  ...
```

## Revision History

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2024-01-01 | Initial version | - |

## References

1. Andrew J. Page, Jacqueline A. Keane. (2018) Rapid multi-locus sequence typing direct from uncorrected long reads using Krocus. PeerJ 6:e5233
2. https://github.com/andrewjpage/krocus
3. https://pubmlst.org/
