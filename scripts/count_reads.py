import argparse
import gzip
import csv
import sys
import time
from collections import Counter

def parse_arguments():
    parser = argparse.ArgumentParser(description='Count sequencing reads mapping to a reference library.')
    parser.add_argument('-r', '--reads', required=True, help='Path to input FASTQ or FASTQ.GZ file')
    parser.add_argument('-l', '--library', required=True, help='Path to reference library CSV (Format: ID, Sequence)')
    parser.add_argument('-o', '--output', required=True, help='Path to output CSV file')
    parser.add_argument('--output-trimmed', help='Optional path to save all trimmed reads and their match status (CSV)')
    parser.add_argument('--trim-start', type=int, default=0, help='Number of bases to trim from start of read (default: 0)')
    parser.add_argument('--trim-length', type=int, default=None, help='Length of sequence to keep after trimming (default: Match library length)')
    parser.add_argument('--revcomp', action='store_true', help='Search for reverse complement of library sequences')
    
    # QC Arguments
    parser.add_argument('--min-base-q', type=int, default=28, help='Minimum Phred quality score for a base to be considered "good" (default: 28)')
    parser.add_argument('--max-bad-freq', type=float, default=10.0, help='Maximum percentage of "bad" bases allowed in a read before dropping it (default: 10%%)')
    
    return parser.parse_args()

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def load_library(filepath):
    """
    Loads the reference library into a dictionary.
    Expects CSV format: ID, Sequence
    """
    library_map = {}
    library_ids = {}
    sequence_lengths = set()
    
    print(f"Loading library from {filepath}...")
    
    try:
        with open(filepath, mode='r') as f:
            reader = csv.reader(f)
            header = next(reader, None) # Skip header if it exists, or handle appropriately
            
            # Heuristic to check if first row is header
            if header and (len(header[1]) < 10 or "seq" in header[1].lower()):
                 # It was likely a header, continue
                 pass
            else:
                # It wasn't a header, process it
                f.seek(0)
                
            for row in reader:
                if len(row) < 2: continue
                seq_id = row[0].strip()
                sequence = row[1].strip().upper()
                
                library_map[sequence] = seq_id
                library_ids[seq_id] = 0 # Initialize count
                sequence_lengths.add(len(sequence))
                
    except Exception as e:
        print(f"Error loading library: {e}")
        sys.exit(1)
        
    print(f"Loaded {len(library_map)} unique sequences.")
    
    # Determine target length for trimming if not specified
    target_len = None
    if len(sequence_lengths) == 1:
        target_len = list(sequence_lengths)[0]
    elif len(sequence_lengths) > 1:
        print(f"Warning: Library contains sequences of varying lengths: {sequence_lengths}")
        print("Exact matching requires reads to be trimmed to the exact library sequence.")
    
    return library_map, library_ids, target_len

def process_reads(reads_path, library_map, library_ids, trim_start, trim_length, revcomp=False, trimmed_output_path=None, min_base_q=28, max_bad_freq=10.0):
    """
    Parses FASTQ and counts matches with QC filtering.
    """
    print(f"Processing reads from {reads_path}...")
    print(f"QC Settings: dropping reads with >={max_bad_freq}% bases under Q{min_base_q}")
    
    # Determine file open mode
    if reads_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt' # read text mode
    else:
        open_func = open
        mode = 'r'

    total_reads = 0
    skipped_qc = 0
    matched_reads = 0
    
    # Pre-calculate QC threshold character
    # Phred+33 is standard for Illumina
    # If ASCII < (min_base_q + 33), the quality is too low.
    qc_char_limit = min_base_q + 33
    max_bad_ratio = max_bad_freq / 100.0
    
    # Open trimmed output file if requested
    trimmed_file = None
    if trimmed_output_path:
        print(f"Saving trimmed reads to {trimmed_output_path}...")
        trimmed_file = open(trimmed_output_path, 'w')
        trimmed_file.write("Read_ID,Trimmed_Sequence,Matched_Library_ID\n")
    
    start_time = time.time()

    try:
        with open_func(reads_path, mode) as f:
            while True:
                # FASTQ format: 4 lines per record
                # Line 1: ID
                identifier = f.readline()
                if not identifier: break # End of file
                
                # Line 2: Sequence
                sequence = f.readline().strip().upper()
                
                # Line 3: +
                f.readline()
                
                # Line 4: Quality
                quality = f.readline().strip()
                
                total_reads += 1
                
                # --- QC STEP ---
                # Count bases with ASCII code < limit
                low_quality_count = 0
                for char in quality:
                    if ord(char) < qc_char_limit:
                        low_quality_count += 1
                
                if len(quality) > 0 and (low_quality_count / len(quality)) >= max_bad_ratio:
                    skipped_qc += 1
                    continue # Skip this read entirely
                # ----------------
                
                # Trim Sequence
                if trim_length:
                    query_seq = sequence[trim_start : trim_start + trim_length]
                else:
                    query_seq = sequence[trim_start:]

                # Apply Reverse Complement if requested (after trimming)
                if revcomp:
                    query_seq = get_reverse_complement(query_seq)
                
                # Check for match
                matched_id = ""
                if query_seq in library_map:
                    seq_id = library_map[query_seq]
                    library_ids[seq_id] += 1
                    matched_reads += 1
                    matched_id = seq_id
                
                # Write to trimmed file if open
                if trimmed_file:
                    # Strip newline and remove '@' from FASTQ header for cleaner CSV
                    clean_id = identifier.strip()[1:]
                    trimmed_file.write(f"{clean_id},{query_seq},{matched_id}\n")
                
                if total_reads % 100000 == 0:
                    print(f"Processed {total_reads} reads... ({matched_reads} matched, {skipped_qc} low qual)", end='\r')
    finally:
        if trimmed_file:
            trimmed_file.close()

    elapsed = time.time() - start_time
    passed_reads = total_reads - skipped_qc
    
    print(f"\nFinished processing {total_reads} reads in {elapsed:.2f} seconds.")
    print(f"Skipped (Low QC): {skipped_qc} ({(skipped_qc/total_reads)*100:.2f}%)")
    print(f"Reads Processed:  {passed_reads}")
    if passed_reads > 0:
        print(f"Total Matched:    {matched_reads} ({(matched_reads/passed_reads)*100:.2f}% of passed reads)")
    
    return library_ids

def write_output(output_path, library_ids):
    print(f"Writing results to {output_path}...")
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Library_ID', 'Count'])
        
        # Sort by count descending
        sorted_results = sorted(library_ids.items(), key=lambda x: x[1], reverse=True)
        
        for seq_id, count in sorted_results:
            writer.writerow([seq_id, count])

def main():
    args = parse_arguments()
    
    library_map, library_ids, lib_seq_len = load_library(args.library)
    
    # Logic to determine cutting length
    cut_len = args.trim_length
    if cut_len is None and lib_seq_len is not None:
        cut_len = lib_seq_len
        print(f"Auto-detected trim length from library: {cut_len}bp")
    
    final_counts = process_reads(
        args.reads, 
        library_map, 
        library_ids, 
        args.trim_start, 
        cut_len, 
        args.revcomp, 
        args.output_trimmed,
        args.min_base_q,
        args.max_bad_freq
    )
    
    write_output(args.output, final_counts)
    print("Done.")

if __name__ == "__main__":
    main()