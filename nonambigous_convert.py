#!/usr/bin/env python3

import sys

"""
Step 2 script for File 1 (A↔C/G or T↔C/G).
Now creates a table with columns:
  SNP_ID, Sample_ID, VCF_REF, VCF_ALT, AlleleA, AlleleB, GenotypeAB, Allele1_AB, Allele2_AB, Strand

Usage:
  python step2_file1.py file1.vcf file1_table.txt
"""

def classify_alleles(ref, alt):
    """
    Determine AlleleA, AlleleB, and strand for a known A↔C/G or T↔C/G site.

    Illumina rules:
      - If {A, C/G}, it's TOP => A=AlleleA, C/G=AlleleB
      - If {T, C/G}, it's BOT => T=AlleleA, C/G=AlleleB
    """
    r = ref.upper()
    a = alt.upper()
    pair = {r, a}

    # A↔C/G => TOP
    # T↔C/G => BOT
    if pair.issubset({'A', 'C', 'G'}):
        # Must be A with C or G => TOP
        strand = "TOP"
        # Allele A = 'A', Allele B = 'C' or 'G'
        if r == 'A':
            alleleA = 'A'
            alleleB = a
        elif a == 'A':
            alleleA = 'A'
            alleleB = r
        else:
            return None, None, None
    else:
        # T↔C/G => BOT
        strand = "BOT"
        # Allele A = 'T', Allele B = 'C' or 'G'
        if r == 'T':
            alleleA = 'T'
            alleleB = a
        elif a == 'T':
            alleleA = 'T'
            alleleB = r
        else:
            return None, None, None

    return alleleA, alleleB, strand

def convert_genotype_to_AB(genotype, ref_base, alt_base, alleleA, alleleB):
    """
    Convert the diploid genotype (e.g. "0|1", "0/1", "1|1") into A/B notation:
      - If index=0 => REF
      - If index=1 => ALT
    Then check whether REF==AlleleA or AlleleB, likewise for ALT.

    Returns:
      (GenotypeAB, Allele1_AB, Allele2_AB)

    Where GenotypeAB is sorted ("AA","AB","BB"), and Allele1_AB, Allele2_AB
    reflect the unphased order in the original call.
    """
    # Missing or invalid genotype
    if '.' in genotype:
        return "NA", "NA", "NA"

    # Normalize phased/unphased calls
    alleles = genotype.replace('|','/').split('/')
    if len(alleles) != 2:
        return "NA", "NA", "NA"

    ab_list = []
    for x in alleles:
        if x == '0':
            # REF => could be AlleleA or AlleleB
            ab_list.append("A" if ref_base == alleleA else "B")
        elif x == '1':
            # ALT => could be AlleleA or AlleleB
            ab_list.append("A" if alt_base == alleleA else "B")
        else:
            return "NA", "NA", "NA"

    # GenotypeAB => sorted, e.g. "AB" for heterozygote
    genotypeAB = "".join(sorted(ab_list))
    allele1_AB, allele2_AB = ab_list[0], ab_list[1]

    return genotypeAB, allele1_AB, allele2_AB

def main():
    if len(sys.argv) < 3:
        print("Usage: python nonambigous_convert.py file1.vcf file1_table.txt", file=sys.stderr)
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_txt = sys.argv[2]

    with open(input_vcf, 'r') as fin, open(output_txt, 'w') as fout:
        # Write header
        fout.write("SNP_ID\tSample_ID\tVCF_REF\tVCF_ALT\tAlleleA\tAlleleB\tGenotypeAB\tAllele1_AB\tAllele2_AB\tStrand\tSourceGT\n")

        sample_names = []
        for line in fin:
            line = line.strip()
            if line.startswith("##"):
                # Skip metadata lines
                continue
            if line.startswith("#CHROM"):
                # Parse sample names from header
                cols = line.split("\t")
                sample_names = cols[9:]
                continue

            cols = line.split("\t")
            if len(cols) < 10:
                # Invalid line
                continue

            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            # Next columns might be QUAL, FILTER, INFO, FORMAT, then genotypes
            format_field = cols[8]
            sample_fields = cols[9:]

            # Build SNP ID
            snp_id = vid if vid != "." else f"{chrom}:{pos}"

            # Classify Allele A/B and Strand
            alleleA, alleleB, strand = classify_alleles(ref, alt)
            if alleleA is None:
                # Not a valid A↔C/G or T↔C/G site
                continue

            # Convert each sample genotype
            for i, sample_data in enumerate(sample_fields):
                sample_id = sample_names[i]
                # Typically "GT:DP:..." => first subfield is GT
                source_gt = sample_data.split(":")[0] if ":" in sample_data else sample_data

                genotypeAB, allele1_AB, allele2_AB = convert_genotype_to_AB(
                    source_gt, ref.upper(), alt.upper(), alleleA, alleleB
                )

                # Write row
                fout.write(
                   f"{snp_id}\t{sample_id}\t{ref}\t{alt}\t"
                   f"{alleleA}\t{alleleB}\t{genotypeAB}\t"
                   f"{allele1_AB}\t{allele2_AB}\t{strand}\t"
                   f"{source_gt}\n"
                   )
                   
if __name__ == "__main__":
    main()
