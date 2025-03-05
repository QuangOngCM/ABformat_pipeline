#!/usr/bin/env python3

import sys
import re

"""
Usage:
  python make_genotype_file_fimpute_dynamic.py \
         genotype_file.txt \
         snp_info_file.txt \
         genotype_file_fimpute.txt

What it does:
1) Reads genotype_file.txt to store sample lines like:
     SampleA 01250
   and to find the maximum length of the string "ID + 1-space + ChipNumber".
2) Reads snp_info_file.txt, e.g.:
     SNP  Chr  Pos  Chip1
   to parse out the chip number (e.g. '1') from the last column name 'Chip1'.
3) Writes genotype_file_fimpute.txt with a header:
     ID Chip Genotypes
   and one line per sample:
     <ID> <space> <ChipNumber> [some spaces] <Genotypes>
   so that all Genotypes start at the same column.

   - Exactly 1 space between ID and Chip
   - Enough spaces after Chip so that the genotype column lines up.
"""

def main():
    if len(sys.argv) < 4:
        print("Usage: python make_genotype_file_fimpute_dynamic.py "
              "genotype_file.txt snp_info_file.txt genotype_file_fimpute.txt",
              file=sys.stderr)
        sys.exit(1)

    genotype_file = sys.argv[1]
    snp_info_file = sys.argv[2]
    output_file   = sys.argv[3]

    # 1) Read genotype_file.txt and store lines.
    #    Each line is "SampleID  <genotypeCodes>".
    #    We'll also find the maximum length of "SampleID + 1space + ChipNumber"
    #    but we don't know ChipNumber yet, so read all sample IDs first.

    lines = []
    sample_ids = []
    with open(genotype_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(" ", 1)
            if len(parts) < 2:
                continue
            sample_id = parts[0]
            geno_str  = parts[1]
            lines.append((sample_id, geno_str))
            sample_ids.append(sample_id)

    # 2) Parse snp_info_file.txt header to extract the chip number from last column (e.g. "Chip1" -> "1").
    with open(snp_info_file, 'r') as f:
        header_line = f.readline().strip()
        cols = header_line.split('\t')  # e.g. ['SNP','Chr','Pos','Chip1']
        last_col = cols[-1]            # e.g. 'Chip1'

    # Extract digits from last_col if it matches something like "Chip1", "Chip2", etc.
    match = re.match(r"Chip(\d+)$", last_col, flags=re.IGNORECASE)
    if match:
        chip_number = match.group(1)  # e.g. '1'
    else:
        chip_number = "1"  # fallback if not matched

    # 3) We still need to figure out how wide "ID + ' ' + chip_number" can be.
    #    Since we didn't know chip_number earlier, let's build that string for each sample now,
    #    then measure its length to find a max.

    max_label_len = 0
    label_list = []
    for (sample_id, geno_str) in lines:
        # "SampleA Chip1" => sample_id + " " + chip_number
        label_str = sample_id + " " + chip_number
        label_len = len(label_str)
        if label_len > max_label_len:
            max_label_len = label_len
        label_list.append((sample_id, geno_str, label_str))

    # 4) Write genotype_file_fimpute.txt
    #    We want a header:
    #      "ID Chip Genotypes"
    #    Then each sample line with:
    #      <sampleID> <space> <chipNumber> <spaces> <genotype_string>
    #    So that <genotype_string> starts in the same column for all lines.

    with open(output_file, 'w') as fout:
        # Write the header
        # The header "ID Chip" itself is 6 chars. We'll find how many spaces we need so that
        # "Genotypes" starts at column (max_label_len + 1).
        # So if max_label_len is 10, we want "Genotypes" to begin at column 11.
        header_left = "ID Chip"  # length=6
        needed_spaces = (max_label_len + 1) - len(header_left)
        if needed_spaces < 1:
            needed_spaces = 1
        header_line = header_left + (" " * needed_spaces) + "Genotypes"
        fout.write(header_line + "\n")

        # Now each sample line
        for (sample_id, geno_str, label_str) in label_list:
            # label_str = sample_id + " " + chip_number
            current_len = len(label_str)
            # We want genotype to start at column (max_label_len + 1)
            spaces_needed = (max_label_len + 1) - current_len
            if spaces_needed < 1:
                spaces_needed = 1
            out_line = label_str + (" " * spaces_needed) + geno_str
            fout.write(out_line + "\n")


if __name__ == "__main__":
    main()

