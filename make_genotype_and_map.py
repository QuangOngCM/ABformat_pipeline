#!/usr/bin/env python3

import sys

"""
Usage:
  python make_genotype_and_map.py final_merged.txt genotype_file.txt snp_map.txt snp_info_file.txt

Outputs:
1) genotype_file.txt
   - One line per sample
   - Format: sampleID <space> genotypeCodesConcatenated
     (AA -> 0, AB -> 1, BB -> 2, missing -> 5)
     in the order SNPs appear in final_merged.txt.

2) snp_map.txt
   - Columns: SNP_ID, CHROM, POS, INDEX
   - INDEX is the 1-based order of first appearance in final_merged.

3) snp_info_file.txt
   - Exactly the same content (rows) as snp_map.txt,
   - but with a different header: "SNP\tChr\tPos\tChip1"
"""

# Convert "AA","AB","BB","NA" => numeric code
GENO_CODE = {
    "AA": 0,
    "AB": 1,
    "BB": 2,
    "NA": 5  # missing/invalid
}

def parse_chrom_pos(snp_id):
    """
    Attempt to parse CHROM and POS from snp_id if in form "chr:pos".
    If not parseable, return (snp_id, "0").
    """
    if ':' in snp_id:
        parts = snp_id.split(':')
        chrom = parts[0]
        pos = parts[1]
        return (chrom, pos)
    else:
        # If it's some other ID (e.g. "rs123"), fallback
        return (snp_id, "0")

def main():
    if len(sys.argv) < 5:
        print("Usage: python make_genotype_and_map.py final_merged.txt genotype_file.txt snp_map.txt snp_info_file.txt", file=sys.stderr)
        sys.exit(1)

    merged_file = sys.argv[1]
    genotype_file = sys.argv[2]
    snp_map_file = sys.argv[3]
    snp_info_file = sys.argv[4]

    # We'll store:
    #   snp_index_list: list of SNP IDs in the order they appear
    #   snp_index_dict: map from snp_id -> index (0-based, we output 1-based)
    snp_index_list = []
    snp_index_dict = {}

    # sample_genotypes[sampleID] = dict of snp_index -> numeric code
    sample_genotypes = {}

    # 1) Read final_merged.txt
    with open(merged_file, 'r') as fin:
        header = fin.readline().rstrip('\n').split('\t')

        # We need columns: SNP_ID, Sample_ID, GenotypeAB
        try:
            snp_id_idx = header.index("SNP_ID")
            sample_id_idx = header.index("Sample_ID")
            genotypeAB_idx = header.index("GenotypeAB")
        except ValueError:
            print("ERROR: final_merged.txt must have columns SNP_ID, Sample_ID, GenotypeAB", file=sys.stderr)
            sys.exit(1)

        for line in fin:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            snp_id = cols[snp_id_idx]
            sample_id = cols[sample_id_idx]
            genoAB = cols[genotypeAB_idx]  # e.g. "AA","AB","BB","NA"

            # Assign index if new SNP
            if snp_id not in snp_index_dict:
                snp_index = len(snp_index_list)
                snp_index_list.append(snp_id)
                snp_index_dict[snp_id] = snp_index
            else:
                snp_index = snp_index_dict[snp_id]

            # Convert genotype to numeric
            geno_code = GENO_CODE.get(genoAB, 5)

            # Store in sample_genotypes
            if sample_id not in sample_genotypes:
                sample_genotypes[sample_id] = {}
            sample_genotypes[sample_id][snp_index] = geno_code

    # 2) Write snp_map.txt
    with open(snp_map_file, 'w') as fout_map:
        fout_map.write("SNP_ID\tCHROM\tPOS\tINDEX\n")
        for i, snp_id in enumerate(snp_index_list):
            chrom, pos = parse_chrom_pos(snp_id)
            index_1based = i + 1
            fout_map.write(f"{snp_id}\t{chrom}\t{pos}\t{index_1based}\n")

    # 3) Write snp_info_file.txt (same rows, different header)
    with open(snp_info_file, 'w') as fout_info:
        fout_info.write("SNP\tChr\tPos\tChip1\n")
        for i, snp_id in enumerate(snp_index_list):
            chrom, pos = parse_chrom_pos(snp_id)
            index_1based = i + 1
            fout_info.write(f"{snp_id}\t{chrom}\t{pos}\t{index_1based}\n")

    # 4) Write genotype_file.txt
    #    - one line per sample => "sampleID <space> genotypeCodes"
    #    - genotypeCodes are concatenated (no spaces between them)
    all_samples = sorted(sample_genotypes.keys())
    with open(genotype_file, 'w') as fout_geno:
        for sample_id in all_samples:
            codes_list = []
            for snp_index in range(len(snp_index_list)):
                code_val = sample_genotypes[sample_id].get(snp_index, 5)
                codes_list.append(str(code_val))
            # Concatenate with no spaces
            genotype_string = "".join(codes_list)
            fout_geno.write(f"{sample_id} {genotype_string}\n")

if __name__ == "__main__":
    main()
