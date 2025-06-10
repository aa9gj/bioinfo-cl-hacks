# Bioinformatics Hacks

This section provides specialized command-line hacks for common bioinformatics tasks, focusing on manipulating and analyzing standard file formats like FASTQ, SAM/BAM, VCF, FASTA, and BED.

## FASTQ Processing

* **Count reads in a FASTQ file (gzipped or not):**
    ```bash
    # For gzipped FASTQ
    zcat input.fastq.gz | wc -l | awk '{print $1/4}'

    # For uncompressed FASTQ
    wc -l input.fastq | awk '{print $1/4}'
    ```
* **Subsample FASTQ files (e.g., 1 million reads):**
    ```bash
    # Using seqtk (recommended for large files)
    # Ensure seqtk is installed and in your PATH, or load via module
    # module load seqtk/1.3

    seqtk sample -s100 read1.fastq.gz 1000000 > subsampled_read1.fastq
    seqtk sample -s100 read2.fastq.gz 1000000 > subsampled_read2.fastq
    ```
* **Basic Quality Trimming with `fastp`:**
    ```bash
    # Ensure fastp is installed and in your PATH, or load via module
    # module load fastp/0.23.2

    fastp -i read1.fastq.gz -o trimmed_read1.fastq.gz \
          -I read2.fastq.gz -O trimmed_read2.fastq.gz \
          --json trimmed_stats.json --html trimmed_stats.html \
          --thread 8
    ```
* **Check FastQ Read Names for consistency (for paired-end reads):**
    ```bash
    # Extract only the read identifiers (e.g., '@readname/1' or '@readname 1:...')
    # and compare sorted lists for R1 and R2 to ensure they match.
    # This is useful for identifying issues before alignment.
    zcat read1.fastq.gz | paste - - - - | cut -f1 | sed 's/ /_/g' | sort > read1_ids.txt
    zcat read2.fastq.gz | paste - - - - | cut -f1 | sed 's/ /_/g' | sort > read2_ids.txt

    diff read1_ids.txt read2_ids.txt

    # If 'diff' produces no output, the IDs match.
    # For speed on large files, consider comparing checksums of the ID files:
    # md5sum read1_ids.txt read2_ids.txt
    ```

---

## SAM/BAM Processing (using `samtools`)

`samtools` is an essential tool for manipulating and analyzing SAM/BAM/CRAM alignment files.

* **Convert SAM to BAM, sort, and index:**
    ```bash
    # Convert SAM to BAM (view -bS)
    samtools view -bS input.sam > output.bam

    # Sort BAM by coordinate (-o specifies output)
    samtools sort output.bam -o output.sorted.bam

    # Index the sorted BAM (creates .bai file)
    samtools index output.sorted.bam

    # In one go (pipe from view to sort, and then index)
    samtools view -bS input.sam | samtools sort -o output.sorted.bam -
    samtools index output.sorted.bam
    ```
* **Count total reads in a BAM file:**
    ```bash
    samtools view -c alignment.bam
    ```
* **Count mapped reads:**
    ```bash
    samtools view -c -F 4 alignment.bam
    # -F 4: Exclude reads where FLAG bit 4 (unmapped) is set
    ```
* **Count unmapped reads:**
    ```bash
    samtools view -c -f 4 alignment.bam
    # -f 4: Include only reads where FLAG bit 4 (unmapped) is set
    ```
* **Extract mapped reads only:**
    ```bash
    samtools view -b -F 4 input.bam > mapped_reads.bam
    ```
* **Extract unmapped reads only:**
    ```bash
    samtools view -b -f 4 input.bam > unmapped_reads.bam
    ```
* **Filter BAM by mapping quality (MAPQ > 20):**
    ```bash
    samtools view -b -q 20 input.bam > high_quality_mapped.bam
    ```
* **Get alignment statistics (`flagstat`):**
    ```bash
    samtools flagstat input.bam
    ```
* **Calculate coverage depth per base (`depth`):**
    ```bash
    samtools depth -a input.sorted.bam > coverage.txt
    # -a: Output all positions (including 0 coverage)
    # Or for specific regions (requires BED file)
    samtools depth -b regions.bed input.sorted.bam > region_coverage.txt
    ```
* **Extract reads from a specific genomic region:**
    ```bash
    samtools view -b input.sorted.bam "chr1:1000-2000" > chr1_region.bam
    ```
* **Merge multiple BAM files:**
    ```bash
    samtools merge merged.bam input1.bam input2.bam input3.bam
    # Requires inputs to be sorted in the same way (usually coordinate)
    ```

---

## VCF Processing (using `bcftools` and `vcftools`)

`bcftools` is a powerful suite for manipulating and analyzing VCF and BCF files. `vcftools` is also widely used for filtering and statistics.

* **View VCF header and first few variants:**
    ```bash
    bcftools view input.vcf.gz | head
    ```
* **Compress and index VCF (essential for many `bcftools` commands):**
    ```bash
    bgzip input.vcf       # Creates input.vcf.gz
    bcftools index input.vcf.gz # Creates input.vcf.gz.csi or .tbi
    ```
* **Filter VCF by variant quality (QUAL > 30):**
    ```bash
    bcftools view -i 'QUAL>30' input.vcf.gz > filtered_qual.vcf
    ```
* **Filter by INFO field depth (DP > 10):**
    ```bash
    bcftools view -i 'INFO/DP>10' input.vcf.gz > filtered_dp.vcf
    ```
* **Filter by FORMAT field genotype quality (GQ >= 20) for all samples:**
    ```bash
    # This filters based on the Genotype Quality for each sample (GQ)
    bcftools filter -i 'FMT/GQ>=20' input.vcf.gz -Oz -o filtered_gq.vcf.gz
    bcftools index filtered_gq.vcf.gz
    ```
* **Extract specific fields from VCF (e.g., CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO/DP, INFO/AF):**
    ```bash
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/DP\t%INFO/AF\n' input.vcf.gz > variants_summary.tsv
    ```
* **Calculate allele frequencies (simple `bcftools` plugin):**
    ```bash
    bcftools +af input.vcf.gz -Oz -o allele_frequencies.vcf.gz
    bcftools index allele_frequencies.vcf.gz
    ```
* **Merge multiple VCF files:**
    ```bash
    bcftools merge file1.vcf.gz file2.vcf.gz file3.vcf.gz -Oz -o merged.vcf.gz
    bcftools index merged.vcf.gz
    ```
* **Annotate VCF with `VEP` (Variant Effect Predictor) - example snippet:**
    ```bash
    # VEP usage is highly specific to your setup (cache, plugins, etc.)
    # This is a very basic example; consult VEP documentation for full details.
    # module load ensembl-vep/109 # Or your installed version

    vep --input input.vcf \
        --output annotated.vcf \
        --cache --dir_cache /path/to/vep_cache \
        --assembly GRCh38 --symbol --tsl --fork 8 \
        --force_overwrite
    ```

---

## FASTA/BED Processing

* **Count sequences in a FASTA file:**
    ```bash
    grep -c ">" input.fasta
    ```
* **Extract specific sequences from a FASTA (using `samtools faidx`):**
    ```bash
    # Index the FASTA first (only needs to be done once per FASTA file)
    samtools faidx reference.fasta

    # Extract a specific region (e.g., chr1 from position 1000 to 2000)
    samtools faidx reference.fasta chr1:1000-2000 > extracted_region.fasta

    # Extract multiple full sequences by name (e.g., from a list of gene IDs)
    echo -e "gene_A\ngene_B\ngene_C" > sequence_names.txt
    xargs samtools faidx reference.fasta < sequence_names.txt > extracted_genes.fasta
    ```
* **Basic BED file manipulation (`bedtools`):**
    `bedtools` is a powerful suite for genomic interval operations.

    ```bash
    # Ensure bedtools is installed and in your PATH, or load via module
    # module load bedtools/2.30.0

    # Intersect two BED files (find overlapping regions)
    bedtools intersect -a file1.bed -b file2.bed > overlapping_regions.bed

    # Merge overlapping/abutting regions in a BED file
    bedtools merge -i input.bed > merged_regions.bed

    # Subtract regions in B from A (remove B from A)
    bedtools subtract -a file_A.bed -b file_B.bed > A_minus_B.bed

    # Find the closest feature in B for each feature in A
    bedtools closest -a query.bed -b target.bed > closest_features.bed

    # Count overlaps (e.g., how many genes overlap with enhancers)
    bedtools intersect -a genes.bed -b enhancers.bed -c > genes_with_enhancer_counts.bed
    ```
