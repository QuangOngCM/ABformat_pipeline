#!/usr/bin/env python3

import sys

"""
make_genotype_files_with_snpinfo.py

Reads:
  1) unified_AB_table_sorted.txt (columns: SNP_ID, Sample_ID, GenotypeAB, Allele1_AB, Allele2_AB)
  2) snp_info_file.txt (columns: SNP_ID, ... , Chip1) -> we only need SNP_ID and the 'Chip1' order

Outputs:
  1) unified_genotype_file.txt:
       One line per sample: "SampleID <space> genotype_string"
       Where genotype_string is built by converting GenotypeAB:
         AA->0, AB->1, BB->2, anything else->5
       The SNP order is the order in snp_info_file.txt (sorted by Chip1).
  
  2) unified_genotype_file_fimpute.txt:
       FImpute-compatible:
         First line: "ID Chip Genotypes"
         Then each sample line, e.g. "SampleA 1    012512..."
         so that the Genotypes column starts at the same position for every sample.
"""

# Mapping from A/B genotype to numeric codes
GENO_MAP = {
    "AA": "0",
    "AB": "1",
    "BB": "2"
    # If not found, default to "5"
}

def main():
    if len(sys.argv) < 5:
        print("Usage: python make_genotype_files.py "
              "unified_AB_table_sorted.txt snp_info_file.txt "
              "unified_genotype_file.txt unified_genotype_file_fimpute.txt",
              file=sys.stderr)
        sys.exit(1)

    unified_file = sys.argv[1]
    snp_info_file = sys.argv[2]
    out_genofile = sys.argv[3]
    out_fimpute = sys.argv[4]

    # 1) Read snp_info_file.txt to get the SNP order
    #    We'll store an ordered list of SNPs (in the order of Chip1).
    #    We assume columns: [SNP_ID, ..., Chip1], or at least "SNP_ID" and "Chip1".
    snp_order = []
    snp_to_index = {}
    with open(snp_info_file, 'r') as f:
        header = f.readline().strip().split()
        try:
            idx_snpid = header.index("SNP")  # or "SNP" if needed
            idx_chip1 = header.index("Chip1")
        except ValueError:
            print("ERROR: snp_info_file.txt must have columns: SNP and Chip1", file=sys.stderr)
            sys.exit(1)

        # We'll read each line, parse Chip1 as an integer, and store (Chip1, SNP_ID)
        tmp_list = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            snp_id = cols[idx_snpid]
            try:
                chip1_val = int(cols[idx_chip1])
            except:
                # skip if not parseable
                continue
            tmp_list.append((chip1_val, snp_id))
        
        # Sort tmp_list by chip1_val
        tmp_list.sort(key=lambda x: x[0])
        # Now build snp_order
        for (chip1_val, snp_id) in tmp_list:
            snp_order.append(snp_id)

    # 2) Read unified_AB_table_sorted.txt to build genotype data
    #    sample_geno[sample_id][snp_id] = numeric_code
    sample_geno = {}
    with open(unified_file, 'r') as f:
        header = f.readline().strip().split("\t")
        try:
            idx_snpid = header.index("SNP_ID")
            idx_sample = header.index("Sample_ID")
            idx_genoAB = header.index("GenotypeAB")
        except ValueError:
            print("ERROR: unified_AB_table_sorted.txt must have columns: "
                  "SNP_ID, Sample_ID, GenotypeAB, Allele1_AB, Allele2_AB",
                  file=sys.stderr)
            sys.exit(1)
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            snp_id = cols[idx_snpid]
            sample_id = cols[idx_sample]
            genoAB = cols[idx_genoAB]
            # Convert genotype
            numeric_code = GENO_MAP.get(genoAB, "5")
            if sample_id not in sample_geno:
                sample_geno[sample_id] = {}
            sample_geno[sample_id][snp_id] = numeric_code

    # 3) Write genotype_file.txt
    #    "SampleID <space> genotype_string"
    with open(out_genofile, 'w') as fout_geno:
        # Sort sample IDs
        all_samples = sorted(sample_geno.keys())
        for sample_id in all_samples:
            codes = []
            for snp_id in snp_order:
                code = sample_geno[sample_id].get(snp_id, "5")
                codes.append(code)
            genotype_str = "".join(codes)
            fout_geno.write(f"{sample_id} {genotype_str}\n")

    # 4) Write genotype_file_fimpute.txt with alignment
    #    Header: "ID Chip Genotypes"
    #    Each sample line: "SampleID 1 <aligned> genotype_str"
    all_samples = sorted(sample_geno.keys())
    # Determine alignment
    max_id_len = max(len(s) for s in all_samples) if all_samples else 0
    k = 3  # spacing factor
    desired_geno_start = max_id_len + 1 + 1 + k  # sampleID + space + "1" + k spaces

    with open(out_fimpute, 'w') as fout_fimp:
        fout_fimp.write("ID Chip Genotypes\n")
        for sample_id in all_samples:
            codes = []
            for snp_id in snp_order:
                code = sample_geno[sample_id].get(snp_id, "5")
                codes.append(code)
            genotype_str = "".join(codes)

            prefix = f"{sample_id} 1"
            spaces_needed = desired_geno_start - len(prefix)
            if spaces_needed < 1:
                spaces_needed = 1
            spaces = " " * spaces_needed
            line_out = prefix + spaces + genotype_str
            fout_fimp.write(line_out + "\n")

    print("Done! Created:")
    print(f"  {out_genofile}")
    print(f"  {out_fimpute}")

if __name__ == "__main__":
    main()
