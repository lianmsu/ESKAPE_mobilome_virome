#!/usr/bin/env python3

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 01_phage_sequence
# Codon usage and sequence analysis of the isolated phages
# Manuscript      : Fig. 6
# Original path   : lfat1t01:70experiments/01phage_sequence/apc.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

"""
Converted from original Perl script:
    Author: Joseph Fass (UC Davis Genome Center)
    This script tests each sequence in a FASTA file for self-overlap using LAST.
    If an overlap is found, it trims the 5' copy, joins ends, and permutes the
    join point to the middle of the output sequence. Optionally aligns reads
    with BWA for validation.
"""

import sys
import os
import re
import subprocess
import argparse
import tempfile
import shutil

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def read_fasta(filename):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:].split()[0]   # store only first word of header
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, ''.join(seq_lines)


def run_last(seq_file, temp_prefix):
    """
    Run LAST self-alignment on the given FASTA file.
    Returns (list_of_lines_from_last, path_to_last_output) or (None, error).
    """
    # Build lastdb
    cmd_db = ['lastdb', temp_prefix, seq_file]
    try:
        subprocess.run(cmd_db, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"LAST database creation failed: {e.stderr.decode()}", file=sys.stderr)
        return None, None

    # Run lastal
    last_out = temp_prefix + '.lastal'
    cmd_aln = [
        'lastal', '-s', '1', '-x', '300', '-f', '0', '-T', '1',
        temp_prefix, seq_file
    ]
    with open(last_out, 'w') as out:
        try:
            subprocess.run(cmd_aln, check=True, stdout=out, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"LAST alignment failed: {e.stderr.decode()}", file=sys.stderr)
            return None, None

    # Parse output: keep only lines that are not comments (not starting with '#')
    lines = []
    with open(last_out, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines.append(line.strip())
    return lines, last_out


def parse_last_overlap(lines):
    """
    Parse the LAST output lines.
    Expects exactly 3 lines. Checks that the first line (full self) is ignored,
    and the second/third are the overlaps. Returns (overlap_length, start1, end1, start2, end2)
    or None if not valid.
    """
    if len(lines) != 3:
        return None

    # Split each line into fields (tab-separated)
    # According to LAST -f 0 format:
    # 0: query name, 1: query length, 2: query start, 3: query end,
    # 4: strand, 5: target name, 6: target length, 7: target start, 8: target end, ...
    # However the original Perl script uses indices:
    #   [0] query name (for equality check)
    #   [1] target name   (for printing)
    #   [2] target start  (0-based) (for printing)
    #   [3] overlap length (for trimming) -> corresponds to query end (if query start==0)
    # We will mimic that behavior exactly.
    fields1 = lines[1].split('\t')
    fields2 = lines[2].split('\t')

    # Make sure we have at least 4 fields
    if len(fields1) < 4 or len(fields2) < 4:
        return None

    # Check query names are equal (self alignment)
    if fields1[0] != fields2[0]:
        return None

    # Extract values (all strings -> convert numbers when needed)
    try:
        # overlap length is taken directly from fields[3] (query end)
        overlap_len = int(fields1[3])
        # For printing: target name, start, end
        target_name1 = fields1[1]
        target_start1 = int(fields1[2]) + 1   # convert to 1-based for output
        target_end1   = target_start1 + overlap_len - 1

        target_name2 = fields2[1]
        target_start2 = int(fields2[2]) + 1
        target_end2   = target_start2 + overlap_len - 1
    except (ValueError, IndexError):
        return None

    return overlap_len, target_name1, target_start1, target_end1, target_name2, target_start2, target_end2


def trim_and_permute(sequence, overlap_len):
    """
    Remove overlap_len bases from the 5' end, then permute the sequence
    so that the join point is moved to the middle.
    Returns the new sequence.
    """
    if overlap_len >= len(sequence):
        # Overlap covers whole sequence – can't trim; return empty or original?
        # In original script, substr would return empty; but that would be odd.
        # We'll return the original sequence unchanged (but this case should rarely occur).
        return sequence
    trimmed = sequence[overlap_len:]       # remove overlapping part at start
    mid = len(trimmed) // 2
    new_seq = trimmed[mid:] + trimmed[:mid]
    return new_seq

def trim_sequence(sequence, overlap_len):
    """
    移除序列5'端的重叠部分（overlap_len个碱基），直接返回修剪后的序列。
    """
    if overlap_len >= len(sequence):
        # 如果重叠长度 ≥ 全长，返回空字符串（但实际不会发生）
        return ""
    trimmed = sequence[overlap_len:]
    return trimmed



def bwa_validation(seq_file, reads_fastq, basename, n):
    """
    Run BWA MEM alignment of reads against the permuted sequence.
    Creates BAM file and index.
    """
    # Base name for output files
    out_prefix = f"{basename}.{n}"
    # Index the permuted sequence
    cmd_index = ['bwa', 'index', seq_file]
    # Redirect stderr to log file
    log_file = f"{out_prefix}.aln.log"
    with open(log_file, 'a') as log:
        subprocess.run(cmd_index, stderr=log)

    # Align reads with BWA MEM
    cmd_mem = ['bwa', 'mem', '-M', seq_file, reads_fastq]
    # Output to SAM, compress with gzip
    sam_gz = f"{out_prefix}.sam.gz"
    with open(sam_gz, 'w') as out_sam_gz, open(log_file, 'a') as log:
        proc_mem = subprocess.Popen(cmd_mem, stdout=subprocess.PIPE, stderr=log)
        # Use gzip to compress
        with subprocess.Popen(['gzip'], stdin=proc_mem.stdout, stdout=out_sam_gz, stderr=log) as gzip_proc:
            proc_mem.stdout.close()
            gzip_proc.communicate()
        proc_mem.wait()

    # Convert SAM to BAM, sort, and index
    # Use samtools view -uhS to convert to BAM, then sort
    bam_sorted = f"{out_prefix}.bam"
    cmd_view = ['samtools', 'view', '-uhS', sam_gz]
    cmd_sort = ['samtools', 'sort', '-o', bam_sorted, '-']
    with open(log_file, 'a') as log:
        # Pipe view into sort
        proc_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=log)
        subprocess.run(cmd_sort, stdin=proc_view.stdout, stderr=log)
        proc_view.wait()

    # Index BAM
    cmd_index_bam = ['samtools', 'index', bam_sorted]
    with open(log_file, 'a') as log:
        subprocess.run(cmd_index_bam, stderr=log)

    # Optionally remove intermediate SAM.gz if disk space is a concern
    # os.remove(sam_gz)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Detect and permute circular sequences (plasmids, bacterial chromosomes) "
                    "by trimming self-overlaps found with LAST.",
        usage="%(prog)s [options] <fasta>"
    )
    parser.add_argument('-b', dest='basename', default='permuted',
                        help="Basename for output files (alphanumeric and underscore only). "
                             "Default: 'permuted'")
    parser.add_argument('-r', dest='reads', metavar='<fastq>',
                        help="High-accuracy long reads (e.g., PacBio corrected) for validation "
                             "with BWA MEM.")
    parser.add_argument('fasta', metavar='<fasta>',
                        help="Input FASTA file containing sequences to test.")
    args = parser.parse_args()

    # Sanity check basename
    if not re.match(r'^\w+$', args.basename):
        print("Error: basename must contain only alphanumeric characters and underscores.",
              file=sys.stderr)
        sys.exit(1)

    # Check that required binaries are in PATH (optional but helpful)
    required = ['lastdb', 'lastal']
    if args.reads:
        required.extend(['bwa', 'samtools'])
    for prog in required:
        if shutil.which(prog) is None:
            print(f"Error: {prog} not found in PATH. Please install and ensure it is accessible.",
                  file=sys.stderr)
            sys.exit(1)

    # Read the FASTA file
    sequences = list(read_fasta(args.fasta))
    if not sequences:
        print("No sequences found in input FASTA.", file=sys.stderr)
        sys.exit(1)

    # Temporary files prefix
    temp_prefix = "temp.apc"

    # Process each sequence
    for idx, (header, sequence) in enumerate(sequences, start=1):
        # Write sequence to a temporary FASTA file
        temp_fasta = f"{temp_prefix}.fa"
        with open(temp_fasta, 'w') as f:
            f.write(f">{header}\n{sequence}\n")

        # Run LAST
        print(f"Processing sequence {idx} ...")
        print("  formatting LAST db ... ", end='', flush=True)
        last_lines, last_out = run_last(temp_fasta, temp_prefix)
        if last_lines is None:
            print("FAILED", file=sys.stderr)
            # Save the LAST output for debugging
            if last_out and os.path.exists(last_out):
                shutil.move(last_out, f"apc_aln_{idx}.txt")
                print(f"  LAST output saved as apc_aln_{idx}.txt", file=sys.stderr)
            continue
        print("running LAST self-alignment ... ", end='', flush=True)

        # Parse the output
        overlap_info = parse_last_overlap(last_lines)
        if overlap_info is None:
            print("irregular LAST results (expected 3 alignments)", file=sys.stderr)
            # Save the output file for inspection
            if os.path.exists(last_out):
                shutil.move(last_out, f"apc_aln_{idx}.txt")
            continue

        overlap_len, tname1, start1, end1, tname2, start2, end2 = overlap_info
        print(f"Overlap detected between {tname1}:{start1}-{end1} and {tname2}:{start2}-{end2}")

        # Trim and permute the sequence
        # new_seq = trim_and_permute(sequence, overlap_len)
        # mid = len(new_seq) // 2

        # Trim the sequence (no permutation)
        new_seq = trim_sequence(sequence, overlap_len)
        if not new_seq:
            print(f"  Warning: overlap length {overlap_len} >= sequence length, output may be empty.", file=sys.stderr)

        # Prepare output file
        out_file = f"{args.basename}.{idx}.fa"
        with open(out_file, 'w') as f:
            f.write(f">{header}\n")
            f.write(new_seq + "\n")
        print(f"  Output written to {out_file}")

        # If reads provided, run BWA validation
        if args.reads:
            print("  indexing permuted sequence ... ", end='', flush=True)
            bwa_validation(out_file, args.reads, args.basename, idx)
            print("  done")

    # Clean up temporary files
    for f in os.listdir('.'):
        if f.startswith(temp_prefix):
            os.remove(f)
    # Also remove any leftover LAST db files (temp.apc.*)
    for ext in ['.prj', '.sds', '.ssp', '.tis']:
        f = temp_prefix + ext
        if os.path.exists(f):
            os.remove(f)


if __name__ == '__main__':
    main()
