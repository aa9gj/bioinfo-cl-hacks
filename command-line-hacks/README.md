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
### Text Processing

#### `grep`: Powerful Pattern Matching

* **Basic search:**
    ```bash
    grep "pattern" file.txt
    ```
* **Case-insensitive search:**
    ```bash
    grep -i "pattern" file.txt
    ```
* **Count occurrences of pattern:**
    ```bash
    grep -c "pattern" file.txt
    ```
* **Show lines NOT matching pattern (invert match):**
    ```bash
    grep -v "pattern" file.txt
    ```
* **Show context (lines before/after/around a match):**
    ```bash
    grep -A 2 "pattern" file.txt # 2 lines after match
    grep -B 3 "pattern" file.txt # 3 lines before match
    grep -C 5 "pattern" file.txt # 5 lines around match
    ```
* **Search recursively in files (without specifying specific directories):**
    ```bash
    grep -r "my_string" . # Searches current directory and subdirectories
    ```
* **Extended regular expressions (ERE) with `grep -E` or `egrep`:**
    ```bash
    grep -E "^[A-Z]+ gene" file.txt # Lines starting with one or more uppercase letters then " gene"
    grep -E "(word1|word2)" file.txt # Lines containing either "word1" or "word2"
    ```
* **Perl-compatible regular expressions (PCRE) with `grep -P`:**
    ```bash
    grep -P "\b[0-9]{3}-\d{2}-\d{4}\b" file.txt # Search for social security number pattern
    ```

#### `awk`: Column-based Processing and Data Manipulation

* **Print specific columns (e.g., 1st and 3rd, default delimiter is whitespace):**
    ```bash
    awk '{print $1, $3}' input.txt
    ```
* **Process CSV files (specify comma delimiter):**
    ```bash
    awk -F',' '{print $1, $5}' input.csv
    ```
* **Filter rows based on a condition (e.g., 3rd column value > 100):**
    ```bash
    awk '$3 > 100 {print $0}' input.txt # Print entire line if condition met
    ```
* **Calculate sum of a column:**
    ```bash
    awk '{s+=$2} END {print s}' input.txt # Sum of values in the 2nd column
    ```
* **Add a header to output and process data:**
    ```bash
    awk 'BEGIN {print "ID\tValue"} {print $1"\t"$2*10}' input.txt
    ```

#### `sed`: Stream Editor (Find and Replace)

* **Replace first occurrence of pattern on each line:**
    ```bash
    sed 's/old_text/new_text/' input.txt
    ```
* **Replace all occurrences of pattern on each line (`g` flag for global):**
    ```bash
    sed 's/old_text/new_text/g' input.txt
    ```
* **Replace in-place (modifies original file, often with a backup):**
    ```bash
    sed -i.bak 's/old_text/new_text/g' input.txt # Creates input.txt.bak as backup
    sed -i 's/old_text/new_text/g' input.txt     # No backup, modifies directly
    ```
* **Delete lines matching a pattern:**
    ```bash
    sed '/pattern_to_delete/d' input.txt
    ```
* **Delete empty lines:**
    ```bash
    sed '/^$/d' input.txt
    ```

#### Other Useful Text Tools

* **`cut`**: Extract sections from each line of files.
    ```bash
    cut -f 2,4 -d ',' file.csv # Extract fields 2 and 4 from comma-separated file
    cut -c 1-5,10- file.txt   # Extract characters from position 1-5 and from 10 to end
    ```
* **`paste`**: Merge lines of files side-by-side.
    ```bash
    paste file1.txt file2.txt > merged_file.txt # Merges line 1 of file1 with line 1 of file2, etc.
    ```
* **`sort`**: Sort lines of text files.
    ```bash
    sort -k3,3n file.txt       # Sort by the 3rd column numerically (-n for numeric)
    sort -r file.txt           # Reverse sort
    sort -u file.txt           # Sort and remove duplicate lines
    ```
* **`uniq`**: Report or omit repeated lines. Works on **consecutive** identical lines.
    ```bash
    sort file.txt | uniq -c    # Sort first, then count unique lines
    ```
* **`head`**: Output the first part of files.
    ```bash
    head -n 10 file.txt        # View the first 10 lines
    ```
* **`tail`**: Output the last part of files.
    ```bash
    tail -n 5 file.txt         # View the last 5 lines
    tail -f access.log         # Follow (tail) a file as it grows (useful for logs)
    ```
* **`wc`**: Print newline, word, and byte counts for files.
    ```bash
    wc -l file.txt             # Count lines in a file
    wc -w file.txt             # Count words
    wc -c file.txt             # Count bytes (characters)
    ```
* **Piping (`|`) and Redirection (`>`, `>>`, `<`, `2>`, `&>`):**
    * `command1 | command2`: Pipe the standard output of `command1` to the standard input of `command2`.
    * `command > file`: Redirect standard output to `file` (overwrites).
    * `command >> file`: Redirect standard output to `file` (appends).
    * `command < file`: Redirect `file` to standard input of `command`.
    * `command 2> error.log`: Redirect standard error to `error.log`.
    * `command &> output_and_error.log`: Redirect both standard output and standard error to one file.
