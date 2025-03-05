#!/usr/bin/env python3

import sys
import pysam

"""
Step 2 for File 2 (A/T or C/G), now with additional columns:
  SNP_ID, Sample_ID, VCF_REF, VCF_ALT, AlleleA, AlleleB, GenotypeAB, Allele1_AB, Allele2_AB, Strand, SurroundingSequence, SourceGT

Process:
1. For each SNP in file2.vcf (A/T or C/G only),
2. Fetch flanking bases from reference.fasta,
3. Determine strand by scanning outward for the first A/T in 5' or 3' direction,
4. Assign Allele A/B accordingly,
5. Convert each sample's genotype to A/B format (returning also the unphased order),
6. Output table with the bracketed surrounding sequence.
"""

def find_strand_for_at_snp(chrom, pos, ref_fasta):
    max_dist = 50
    seq_len = ref_fasta.get_reference_length(chrom)
    for dist in range(1, max_dist + 1):
        up_coord = pos - dist
        down_coord = pos + dist
        if up_coord > 0:
            upstream_base = ref_fasta.fetch(chrom, up_coord - 1, up_coord).upper()
            if upstream_base in ['A', 'T']:
                return "TOP"
        if down_coord <= seq_len:
            downstream_base = ref_fasta.fetch(chrom, down_coord - 1, down_coord).upper()
            if downstream_base in ['A', 'T']:
                return "BOT"
    return None

def find_strand_for_cg_snp(chrom, pos, ref_fasta):
    # Same scanning logic applies.
    return find_strand_for_at_snp(chrom, pos, ref_fasta)

def assign_alleles_for_ambiguous(ref_base, alt_base, strand):
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    baseset = {ref_base, alt_base}
    if baseset == {'A', 'T'}:
        return ('A','T') if strand == "TOP" else ('T','A')
    elif baseset == {'C', 'G'}:
        return ('C','G') if strand == "TOP" else ('G','C')
    else:
        return (None, None)

def convert_genotype_to_AB(gt_str, ref_base, alt_base, alleleA, alleleB):
    """
    Converts a genotype string (e.g. "0|1", "0/0", etc.) into A/B format.
    Returns a triple: (GenotypeAB, Allele1_AB, Allele2_AB)
    where GenotypeAB is sorted (e.g. "AB") and Allele1_AB and Allele2_AB reflect the original allele order.
    """
    if '.' in gt_str:
        return "NA", "NA", "NA"
    parts = gt_str.replace('|','/').split('/')
    if len(parts) != 2:
        return "NA", "NA", "NA"
    ab_list = []
    for x in parts:
        if x == '0':
            # If ref equals AlleleA, then allele call is "A", else "B"
            ab_list.append('A' if ref_base.upper() == alleleA else 'B')
        elif x == '1':
            ab_list.append('A' if alt_base.upper() == alleleA else 'B')
        else:
            return "NA", "NA", "NA"
    genotypeAB = "".join(sorted(ab_list))
    allele1_AB = ab_list[0]
    allele2_AB = ab_list[1]
    return genotypeAB, allele1_AB, allele2_AB

def highlight_snp_in_context(context_seq, snp_relpos, allele_set):
    """
    Inserts bracket notation [X/Y] into the context sequence at the SNP position.
    For example, if allele_set is {'C','T'}, it returns something like: ...GCA[C/T]GTA...
    """
    if len(allele_set) == 2:
        label = "/".join(sorted(allele_set))
    else:
        label = "/".join(allele_set)
    return context_seq[:snp_relpos] + "[" + label + "]" + context_seq[snp_relpos+1:]

def main():
    if len(sys.argv) < 4:
        print("Usage: python ambiguous_convert.py file2.vcf reference.fasta file2_table.txt", file=sys.stderr)
        sys.exit(1)
    vcf_file = sys.argv[1]
    ref_fasta_path = sys.argv[2]
    output_file = sys.argv[3]
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    with open(vcf_file, 'r') as fin, open(output_file, 'w') as fout:
        # Write header with all columns
        fout.write("SNP_ID\tSample_ID\tVCF_REF\tVCF_ALT\tAlleleA\tAlleleB\tGenotypeAB\tAllele1_AB\tAllele2_AB\tStrand\tSurroundingSequence\tSourceGT\n")
        sample_names = []
        for line in fin:
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_cols = line.split("\t")
                sample_names = header_cols[9:]
                continue
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            chrom, pos_str, vid, ref_base, alt_base = cols[0], cols[1], cols[2], cols[3], cols[4]
            pos = int(pos_str)
            sample_fields = cols[9:]
            snp_id = vid if vid != "." else f"{chrom}:{pos}"
            # Determine strand using flanking sequences
            base_set = {ref_base.upper(), alt_base.upper()}
            if base_set == {'A','T'}:
                strand = find_strand_for_at_snp(chrom, pos, ref_fasta)
            elif base_set == {'C','G'}:
                strand = find_strand_for_cg_snp(chrom, pos, ref_fasta)
            else:
                continue
            if strand is None:
                continue
            alleleA, alleleB = assign_alleles_for_ambiguous(ref_base, alt_base, strand)
            if alleleA is None:
                continue
            # Get surrounding context (±10 bp)
            flank_size = 10
            chrom_len = ref_fasta.get_reference_length(chrom)
            fetch_start = max(0, pos - 1 - flank_size)
            fetch_end = min(chrom_len, pos - 1 + flank_size + 1)
            context_seq = ref_fasta.fetch(chrom, fetch_start, fetch_end).upper()
            snp_relpos = (pos - 1) - fetch_start
            bracketed_seq = highlight_snp_in_context(context_seq, snp_relpos, base_set)
            for i, sample_data in enumerate(sample_fields):
                sample_id = sample_names[i]
                subfields = sample_data.split(":")
                source_gt = subfields[0]
                genotypeAB, allele1_AB, allele2_AB = convert_genotype_to_AB(source_gt, ref_base, alt_base, alleleA, alleleB)
                fout.write(f"{snp_id}\t{sample_id}\t{ref_base}\t{alt_base}\t{alleleA}\t{alleleB}\t{genotypeAB}\t{allele1_AB}\t{allele2_AB}\t{strand}\t{bracketed_seq}\t{source_gt}\n")
    ref_fasta.close()

if __name__ == "__main__":
    main()
