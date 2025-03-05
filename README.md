    AB Format Pipeline
    This repository contains a suite of Python and Bash scripts that convert raw SNP data in VCF format into Illumina A/B notation, merge ambiguous and non-ambiguous SNPs, and produce final genotype and mapping files. The pipeline includes:

        Split VCF: Separates non-ambiguous (A↔C/G, T↔C/G) from ambiguous (A↔T, C↔G) SNPs.
        File 1 & File 2 Converters: Converts each subset into Illumina A/B format, including strand assignment (TOP/BOT) and source genotype tracking.
        Merge: Combines the outputs into a single table (final_merged.txt).
        Genotype & Map Creation: Generates numeric genotype files (0,1,2,5) and SNP map files for downstream analyses or imputation.

    Key Scripts:

        split_vcf.py: Splits test.vcf into file1.vcf (non-ambiguous) and file2.vcf (ambiguous).
        nonambigous_convert.py: Converts non-ambiguous SNPs to A/B format.
        ambiguous_convert.py: Handles ambiguous SNPs, determining strand from flanking sequences.
        merge_ab_tables.py: Merges tables into a single output.
        make_genotype_and_map.py: Creates numeric genotype and SNP map files from final_merged.txt.
        VCF2ABformat.sh: A one-step script that orchestrates the entire workflow.

    Usage:

        Clone this repository
        Provide a reference FASTA and VCF file
        Run VCF2ABformat.sh <input_vcf> <reference_fasta> #the *.fai index file should be at the same directory with the fasta file
        Check final outputs: final_merged.txt, genotype_file.txt, snp_map.txt, etc.
