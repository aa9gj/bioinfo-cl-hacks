# SLURM Hacks

This section provides a collection of useful SLURM job submission scripts, templates, and tips to help you manage your bioinformatics workloads on HPC clusters.

## Basic SLURM Commands

* `sbatch <script.sh>`: Submit a job script to the SLURM scheduler.
* `squeue`: View currently running and pending jobs.
* `scancel <job_id>`: Cancel a submitted job.
* `sinfo`: Display information about partitions and nodes.
* `sacct -j <job_id>`: Get detailed accounting information for a completed job.

## Key SLURM Directives (`#SBATCH`)

These directives are placed at the top of your job script and tell SLURM how to run your job.

* `#SBATCH --job-name=<name>`: A meaningful name for your job.
* `#SBATCH --nodes=<N>`: Number of nodes to request (usually 1 for most bioinformatics tasks).
* `#SBATCH --ntasks-per-node=<N>`: Number of tasks (processes) to run per node.
* `#SBATCH --cpus-per-task=<N>`: Number of CPU cores requested per task.
* `#SBATCH --mem=<amount>`: Memory requested (e.g., `16G` for 16 Gigabytes).
* `#SBATCH --time=<HH:MM:SS>`: Wall-clock time limit for the job.
* `#SBATCH --output=slurm-%j.out`: Path to the standard output file. `%j` is replaced by the job ID.
* `#SBATCH --error=slurm-%j.err`: Path to the standard error file.
* `#SBATCH --partition=<partition_name>`: Specify a SLURM partition (queue).

---

### SLURM Array Jobs

Array jobs are incredibly useful for running many similar, independent tasks in parallel.

#### `slurm_array_job_template.sh`

A comprehensive template for processing multiple files or samples using `$SLURM_ARRAY_TASK_ID`.

```bash
#!/bin/bash
#SBATCH --job-name=array_process     # Job name
#SBATCH --output=array_%A_%a.out     # Standard output and error log
#SBATCH --error=array_%A_%a.err      # %A is job ID, %a is array task ID
#SBATCH --time=01:00:00              # Wall clock time limit
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=8G                     # Memory per task
#SBATCH --array=1-10                 # Array range (e.g., 1 to 10 tasks)
#SBATCH --export=ALL                 # Export all current environment variables

# Load necessary modules
module load your_bio_software/version

# Define a list of input files or parameters
# This list can be generated dynamically or stored in a file
declare -a SAMPLES=(
    "sample_A"
    "sample_B"
    "sample_C"
    "sample_D"
    "sample_E"
    "sample_F"
    "sample_G"
    "sample_H"
    "sample_I"
    "sample_J"
)

# Get the sample name for the current array task
# SLURM_ARRAY_TASK_ID is 1-indexed, so adjust for 0-indexed array
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting task $SLURM_ARRAY_TASK_ID for sample: $SAMPLE_NAME"
echo "Running on host: $(hostname)"
echo "Current directory: $(pwd)"

# --- Your processing command here ---
# Example: Run a hypothetical analysis on each sample
your_analysis_tool --input_file "${SAMPLE_NAME}.fastq.gz" \
                   --output_prefix "${SAMPLE_NAME}_processed" \
                   --threads $SLURM_CPUS_PER_TASK \
                   --memory $(echo "$SLURM_MEM" | sed 's/G//') # Convert G to integer if needed

if [ $? -eq 0 ]; then
    echo "Task $SLURM_ARRAY_TASK_ID for $SAMPLE_NAME completed successfully."
else
    echo "Task $SLURM_ARRAY_TASK_ID for $SAMPLE_NAME failed."
    exit 1
fi
```
### Using an Input File for Array Jobs

For very long lists of samples, it's better to store them in a file. This allows for dynamic array sizing based on the number of lines in your input list.

#### `slurm_array_job_from_file.sh`

```bash
#!/bin/bash
#SBATCH --job-name=array_file_process
#SBATCH --output=array_file_%A_%a.out
#SBATCH --error=array_file_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-$(wc -l < "samples.txt") # Dynamically set array size based on lines in file
#SBATCH --export=ALL

# Define the file containing your list of samples/inputs, one per line
SAMPLE_LIST_FILE="samples.txt"

# Check if the sample list file exists
if [ ! -f "$SAMPLE_LIST_FILE" ]; then
    echo "Error: Sample list file '<span class="math-inline">SAMPLE\_LIST\_FILE' not found\!"
exit 1
fi
\# Read the specific line corresponding to the current array task
\# sed \-n "</span>{SLURM_ARRAY_TASK_ID}p" reads the Nth line of the file
SAMPLE_IDENTIFIER=<span class="math-inline">\(sed \-n "</span>{SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST_FILE")

# --- Basic logging ---
echo "--- SLURM Array Task ID: $SLURM_ARRAY_TASK_ID ---"
echo "Processing identifier: $SAMPLE_IDENTIFIER"
echo "Start time: <span class="math-inline">\(date\)"
\# \-\-\- Load modules \-\-\-
module load bwa/0\.7\.17
module load samtools/1\.16\.1
\# \-\-\- Your processing command here \-\-\-
\# Example\: Aligning a sample where SAMPLE\_IDENTIFIER is a FASTQ prefix
\# Assumes input files are like 'data/sample\_id\_R1\.fastq\.gz' and 'data/sample\_id\_R2\.fastq\.gz'
READ1\="data/</span>{SAMPLE_IDENTIFIER}_R1.fastq.gz"
READ2="data/<span class="math-inline">\{SAMPLE\_IDENTIFIER\}\_R2\.fastq\.gz"
REFERENCE\="/path/to/genome/reference\.fasta"
OUTPUT\_BAM\="</span>{SAMPLE_IDENTIFIER}.sorted.bam"

echo "Aligning $READ1 and $READ2..."
bwa mem -t $SLURM_CPUS_PER_TASK $REFERENCE $READ1 <span class="math-inline">READ2 \| \\
samtools view \-Sb \- \> "</span>{SAMPLE_IDENTIFIER}.bam"

# Sort and index the BAM file
samtools sort "${SAMPLE_IDENTIFIER}.bam" -o "$OUTPUT_BAM"
samtools index "<span class="math-inline">OUTPUT\_BAM"
\# Clean up unsorted BAM to save space
rm "</span>{SAMPLE_IDENTIFIER}.bam"

if [ $? -eq 0 ]; then
    echo "Task $SLURM_ARRAY_TASK_ID ($SAMPLE_IDENTIFIER) completed successfully."
else
    echo "ERROR: Task $SLURM_ARRAY_TASK_ID ($SAMPLE_IDENTIFIER) failed."
    exit 1
fi
```
---

## Job Dependencies and Chaining

Chain jobs to ensure a workflow progresses sequentially, with each step starting only after the previous one successfully completes.

### Chaining Jobs with `--dependency=afterok` (`chaining_slurm_jobs.sh`)

```bash
#!/bin/bash
# This script submits a series of dependent jobs.
# Submit this script itself to SLURM: sbatch chaining_slurm_jobs.sh

# --- Part 1: Alignment Job ---
echo "Submitting Part 1: Alignment Job..."
# The --parsable option ensures only the Job ID is printed to stdout.
JOB_ID_1=<span class="math-inline">\(sbatch \-\-parsable <<EOF
\#\!/bin/bash
\#SBATCH \-\-job\-name\=align\_step
\#SBATCH \-\-output\=align\_%j\.out
\#SBATCH \-\-error\=align\_%j\.err
\#SBATCH \-\-time\=01\:00\:00
\#SBATCH \-\-cpus\-per\-task\=8
\#SBATCH \-\-mem\=16G
echo "Starting alignment for SampleA\.\.\."
module load bwa/0\.7\.17 samtools/1\.16\.1
bwa mem \-t \\$SLURM\_CPUS\_PER\_TASK /path/to/reference\.fasta /path/to/SampleA\_R1\.fastq\.gz /path/to/SampleA\_R2\.fastq\.gz \| \\
samtools view \-Sb \- \> SampleA\.bam
samtools sort SampleA\.bam \-o SampleA\.sorted\.bam
samtools index SampleA\.sorted\.bam
rm SampleA\.bam \# Clean up unsorted BAM
if \[ \\</span>? -eq 0 ]; then
    echo "Alignment for SampleA completed successfully."
else
    echo "Alignment for SampleA failed."
    exit 1
fi
EOF
)

if [ -z "$JOB_ID_1" ]; then
    echo "ERROR: Failed to submit Part 1 job."
    exit 1
fi
echo "Part 1 (Alignment) submitted with ID: $JOB_ID_1"

# --- Part 2: Variant Calling Job (depends on Part 1 success) ---
echo "Submitting Part 2: Variant Calling Job (dependent on <span class="math-inline">JOB\_ID\_1\)\.\.\."
JOB\_ID\_2\=</span>(sbatch --parsable --dependency=afterok:$JOB_ID_1 <<EOF
#!/bin/bash
#SBATCH --job-name=variant_call_step
#SBATCH --output=variant_call_%j.out
#SBATCH --error=variant_call_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

echo "Starting variant calling for SampleA (dependent on alignment job <span class="math-inline">JOB\_ID\_1\)\.\.\."
module load bcftools/1\.16
samtools mpileup \-Ou \-f /path/to/reference\.fasta SampleA\.sorted\.bam \| \\
bcftools call \-mv \-Oz \-o SampleA\.vcf\.gz
bcftools index SampleA\.vcf\.gz
if \[ \\</span>? -eq 0 ]; then
    echo "Variant calling for SampleA completed successfully."
else
    echo "Variant calling for SampleA failed."
    exit 1
fi
EOF
)

if [ -z "$JOB_ID_2" ]; then
    echo "ERROR: Failed to submit Part 2 job."
    exit 1
fi
echo "Part 2 (Variant Calling) submitted with ID: $JOB_ID_2"

# --- Part 3: Annotation Job (depends on Part 2 success) ---
echo "Submitting Part 3: Annotation Job (dependent on <span class="math-inline">JOB\_ID\_2\)\.\.\."
JOB\_ID\_3\=</span>(sbatch --parsable --dependency=afterok:$JOB_ID_2 <<EOF
#!/bin/bash
#SBATCH --job-name=annotation_step
#SBATCH --output=annotation_%j.out
#SBATCH --error=annotation_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

echo "Starting annotation for SampleA (dependent on variant calling job <span class="math-inline">JOB\_ID\_2\)\.\.\."
module load ensembl\-vep/109 \# Example VEP module, adjust as needed
\# This is a simplified VEP command\. Real VEP usage is more complex\.
\# Consult VEP documentation for full details and cache setup\.
\# vep \-\-input SampleA\.vcf\.gz \-\-output SampleA\.annotated\.vcf \-\-cache \-\-fork \\$SLURM\_CPUS\_PER\_TASK
echo "Annotation for SampleA completed successfully \(simulated\)\." \# Replace with actual command
if \[ \\</span>? -eq 0 ]; then
    echo "Annotation for SampleA completed successfully."
else
    echo "Annotation for SampleA failed."
    exit 1
fi
EOF
)

if [ -z "$JOB_ID_3" ]; then
    echo "ERROR: Failed to submit Part 3 job."
    exit 1
fi
echo "Part 3 (Annotation) submitted with ID: $JOB_ID_3"

echo "All workflow jobs submitted. Monitor with 'squeue' or 'sacct'."
```
