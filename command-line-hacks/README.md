# Command Line Hacks

A collection of general command-line snippets, scripts, and tips that can significantly boost your productivity and mastery of the command line, applicable across various domains.

## Useful Aliases

Add these to your `~/.bashrc` (or `~/.zshrc` or similar shell configuration file) and then `source ~/.bashrc` to make them active in your current session.

* `alias ll='ls -lhtr'`: List files in long format, human-readable sizes, sorted by modification time (reverse).
* `alias duf='du -sh * | sort -rh'`: Show disk usage of current directory contents, sorted by size (human-readable, reverse).
* `alias grep_r='grep -r --exclude-dir={.git,node_modules,build}'`: Recursive grep ignoring common development/build directories.
* `alias untar='tar -xvzf'`: Quick untar command for `.tar.gz` files.
* `alias myip='curl ifconfig.me'`: Get your external IP address from the command line.
* `alias c='clear'`: Clear the terminal screen.
* `alias hg='history | grep '`: Quickly search your command history.

---

## Bash Scripting Tips

### Basic Control Flow

```bash
# If-else statement
# Checks if a file exists
if [ -f "my_file.txt" ]; then
    echo "File exists."
else
    echo "File does not exist."
fi

# For loop through a hardcoded list of items
for item in item1 item2 item3; do
    echo "Processing $item"
done

# For loop through files matching a pattern
# Processes all .fastq.gz files in the current directory
for f in *.fastq.gz; do
    echo "Processing file: $f"
done

# While loop, reading line by line from a file
# Reads a file named 'list.txt' and processes each line
while IFS= read -r line; do
    echo "Processing line: $line"
done < "list.txt"
```
Markdown

# Command Line Hacks

A collection of general command-line snippets, scripts, and tips that can significantly boost your productivity and mastery of the command line, applicable across various domains.

## Useful Aliases

Add these to your `~/.bashrc` (or `~/.zshrc` or similar shell configuration file) and then `source ~/.bashrc` to make them active in your current session.

* `alias ll='ls -lhtr'`: List files in long format, human-readable sizes, sorted by modification time (reverse).
* `alias duf='du -sh * | sort -rh'`: Show disk usage of current directory contents, sorted by size (human-readable, reverse).
* `alias grep_r='grep -r --exclude-dir={.git,node_modules,build}'`: Recursive grep ignoring common development/build directories.
* `alias untar='tar -xvzf'`: Quick untar command for `.tar.gz` files.
* `alias myip='curl ifconfig.me'`: Get your external IP address from the command line.
* `alias c='clear'`: Clear the terminal screen.
* `alias hg='history | grep '`: Quickly search your command history.

---

## Bash Scripting Tips

### Basic Control Flow

```bash
# If-else statement
# Checks if a file exists
if [ -f "my_file.txt" ]; then
    echo "File exists."
else
    echo "File does not exist."
fi

# For loop through a hardcoded list of items
for item in item1 item2 item3; do
    echo "Processing $item"
done

# For loop through files matching a pattern
# Processes all .fastq.gz files in the current directory
for f in *.fastq.gz; do
    echo "Processing file: $f"
done

# While loop, reading line by line from a file
# Reads a file named 'list.txt' and processes each line
while IFS= read -r line; do
    echo "Processing line: $line"
done < "list.txt"

### Error Handling
Include these lines at the top of your scripts for robust error handling. They make your scripts fail fast and explicitly, which is crucial for debugging.
```bash
#!/bin/bash
set -euo pipefail

# -e: Exit immediately if a command exits with a non-zero status.
# -u: Treat unset variables as an error when substituting.
# -o pipefail: The return value of a pipeline is the status of the last command
#              to exit with a non-zero status, rather than the last command in the pipeline.
```
### File Manipulation

#### Finding Files

* **Find files by name:**
    ```bash
    find . -name "*.txt"                    # Find all .txt files in current dir and subdirs
    find /path/to/search -iname "*report*.csv" # Case-insensitive search for CSV files containing "report"
    ```
* **Find files by type (f=file, d=directory, l=symlink):**
    ```bash
    find . -type f -name "*.sh"             # Find all shell scripts
    find . -type d -name "results_*"        # Find directories starting with "results_"
    ```
* **Find files by size:**
    ```bash
    find . -size +1G                        # Files larger than 1GB
    find . -size -10M                       # Files smaller than 10MB
    find . -size +100k -size -1M            # Files between 100KB and 1MB
    ```
* **Find and execute command on found files:**
    ```bash
    find . -name "*.log" -exec rm {} \;      # Delete all .log files (less efficient for many files)
    find . -name "*.tmp" -print0 | xargs -0 rm # Delete .tmp files efficiently (handles spaces in names)
    ```

#### Archiving and Compression

* **Compress a directory into a `.tar.gz` archive:**
    ```bash
    tar -czvf my_archive.tar.gz my_directory/
    ```
* **Decompress a `.tar.gz` file:**
    ```bash
    tar -xzvf my_archive.tar.gz
    ```
* **Compress a single file with gzip:**
    ```bash
    gzip my_file.txt # Creates my_file.txt.gz
    gunzip my_file.txt.gz # Decompresses it back
    ```
* **Compress a single file with bzip2 (often better compression than gzip):**
    ```bash
    bzip2 my_file.txt # Creates my_file.txt.bz2
    bunzip2 my_file.txt.bz2 # Decompresses it back
    ```
