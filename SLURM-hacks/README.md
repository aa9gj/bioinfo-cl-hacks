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
