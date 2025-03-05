#!/usr/bin/env python3

import sys

"""
Merge script ensuring SurroudingSequence is column #11, SourceGT is column #12
"""

def merge_files(file1, file2, output):
    # Final header with 12 columns, where col #11=SurroundingSequence, col #12=SourceGT
    header = (
        "SNP_ID\tSample_ID\tVCF_REF\tVCF_ALT\t"
        "AlleleA\tAlleleB\tGenotypeAB\tAllele1_AB\tAllele2_AB\tStrand\t"
        "SurroundingSequence\tSourceGT\n"
    )

    with open(output, 'w') as fout:
        fout.write(header)

        # Process File1 (which has no SurroudingSequence but does have SourceGT as last col)
        with open(file1, 'r') as f1:
            f1.readline()  # skip file1 header
            for line in f1:
                line = line.strip()
                if not line:
                    continue
                # File1 has 11 columns => the last one is SourceGT
                cols = line.split('\t')
                if len(cols) < 11:
                    continue
                # Remove the last col => source_gt
                source_gt = cols.pop()  # now we have 10 columns
                # Insert "NA" for SurroudingSequence as the second-last column
                cols.append("NA")       # col #11 => SurroudingSequence
                cols.append(source_gt)  # col #12 => SourceGT
                fout.write("\t".join(cols) + "\n")

        # Process File2 (which has all 12 columns including SurroudingSequence & SourceGT)
        with open(file2, 'r') as f2:
            f2.readline()  # skip file2 header
            for line in f2:
                fout.write(line)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python merge_ab_tables.py file1_table.txt file2_table.txt final_merged.txt", file=sys.stderr)
        sys.exit(1)

    merge_files(sys.argv[1], sys.argv[2], sys.argv[3])
