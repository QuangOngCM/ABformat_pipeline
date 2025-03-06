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
   - Exactly the same rows as snp_map.txt,
   - but with a different header: "SNP\tChr\tPos\tChip1"

Logic changes:
- If SNP_ID has the form "chr:pos", we parse from that.
- Otherwise, we retrieve CHR, POS from the columns in final_merged.txt.
"""

# Convert "AA","AB","BB","NA" => numeric code
GENO_CODE = {
    "AA": 0,
    "AB": 1,
    "BB": 2,
    "NA": 5  # missing/invalid
}

def parse_snp_info(snp_id, chr_col, pos_col):
    """
    If snp_id has a colon, parse from snp_id -> (chrom, pos).
    Otherwise, use the columns 'chr_col', 'pos_col'.
    """
    if ':' in snp_id:
        # e.g. "chr1:12345"
        parts = snp_id.split(':', 1)
        chrom = parts[0]
        pos = parts[1]
        return (chrom, pos)
    else:
        # fallback: use the CHR, POS from the file
        return (chr_col, pos_col)

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

    # We'll also store snp_meta[snp_id] = (chrom, pos) for each SNP
    snp_meta = {}

    # 1) Read final_merged.txt
    with open(merged_file, 'r') as fin:
        header = fin.readline().rstrip('\n').split('\t')

        # We need columns: SNP_ID, CHR, POS, Sample_ID, GenotypeAB
        try:
            idx_snp = header.index("SNP_ID")
            idx_chr = header.index("CHR")
            idx_pos = header.index("POS")
            idx_sample = header.index("Sample_ID")
            idx_genoAB = header.index("GenotypeAB")
        except ValueError:
            print("ERROR: final_merged.txt must have columns: SNP_ID, CHR, POS, Sample_ID, GenotypeAB", file=sys.stderr)
            sys.exit(1)

        for line in fin:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                continue

            snp_id = cols[idx_snp]
            chr_val = cols[idx_chr]
            pos_val = cols[idx_pos]
            sample_id = cols[idx_sample]
            genoAB = cols[idx_genoAB]  # e.g. "AA","AB","BB","NA"

            # If new SNP, assign index
            if snp_id not in snp_index_dict:
                snp_index = len(snp_index_list)
                snp_index_list.append(snp_id)
                snp_index_dict[snp_id] = snp_index

                # Determine final chrom, pos from either snp_id or columns
                chrom_parsed, pos_parsed = parse_snp_info(snp_id, chr_val, pos_val)
                snp_meta[snp_id] = (chrom_parsed, pos_parsed)
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
            index_1based = i + 1
            chrom, pos = snp_meta[snp_id]
            fout_map.write(f"{snp_id}\t{chrom}\t{pos}\t{index_1based}\n")

    # 3) Write snp_info_file.txt (same rows, different header)
    with open(snp_info_file, 'w') as fout_info:
        fout_info.write("SNP\tChr\tPos\tChip1\n")
        for i, snp_id in enumerate(snp_index_list):
            index_1based = i + 1
            chrom, pos = snp_meta[snp_id]
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
