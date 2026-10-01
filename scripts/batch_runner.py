import os
import glob
import argparse
import subprocess
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description='Batch process Illumina reads using count_reads.py')
    
    # path args
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing FASTQ files')
    parser.add_argument('-o', '--output-dir', required=True, help='Directory to save output CSVs')
    parser.add_argument('-s', '--script', default='count_reads.py', help='Path to the count_reads.py script (default: ./count_reads.py)')
    
    # arguments to pass through to the main script
    parser.add_argument('-l', '--library', required=True, help='Path to reference library CSV')
    parser.add_argument('--trim-start', type=int, default=0, help='Bases to trim from start')
    parser.add_argument('--trim-length', type=int, default=None, help='Length to keep (optional)')
    parser.add_argument('--revcomp', action='store_true', help='Search for reverse complement')
    parser.add_argument('--min-base-q', type=int, default=28, help='Min Base Quality (QC)')
    parser.add_argument('--max-bad-freq', type=float, default=10.0, help='Max Bad Frequency (QC)')
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # 1. Verify the worker script exists
    if not os.path.isfile(args.script):
        print(f"Error: Could not find the script '{args.script}'")
        sys.exit(1)

    # 2. Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        print(f"Creating output directory: {args.output_dir}")
        os.makedirs(args.output_dir)

    # 3. Find all R1_001.fastq.gz files
    # Search pattern: /path/to/input/*R1_001.fastq.gz
    search_pattern = os.path.join(args.input_dir, "*R1_001.fastq.gz")
    files = glob.glob(search_pattern)
    
    if not files:
        print(f"No files matching '*R1_001.fastq.gz' found in {args.input_dir}")
        sys.exit(1)
        
    print(f"Found {len(files)} samples to process.\n")

    # 4. Iterate and Execute
    for filepath in sorted(files):
        filename = os.path.basename(filepath)
        
        # Extract sample name (remove .fastq.gz and _R1_001 if desired)
        # Strategy: Split by _R1_001 to get the prefix
        sample_name = filename.split('_R1_001')[0]
        
        output_csv = os.path.join(args.output_dir, f"{sample_name}_counts.csv")
        trimmed_csv = os.path.join(args.output_dir, f"{sample_name}_trimmed.csv")
        
        print(f"--> Processing {sample_name}...")
        
        # Construct command
        cmd = [
            sys.executable, args.script,
            '-r', filepath,
            '-l', args.library,
            '-o', output_csv,
            '--output-trimmed', trimmed_csv, # We automatically save trimmed logs for batch runs
            '--trim-start', str(args.trim_start),
            '--min-base-q', str(args.min_base_q),
            '--max-bad-freq', str(args.max_bad_freq)
        ]
        
        # Add optional flags
        if args.trim_length:
            cmd.extend(['--trim-length', str(args.trim_length)])
            
        if args.revcomp:
            cmd.append('--revcomp')

        # Run the command
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {filename}: {e}")
            # We choose not to exit here so one bad file doesn't stop the whole batch
            continue
            
    print("\nBatch processing complete.")

if __name__ == "__main__":
    main()
    