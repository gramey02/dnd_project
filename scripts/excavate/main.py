#!/usr/bin/env python3

"""
This script runs EXCAVATE, a pipeline that generates guide RNA libraries for a given locus, with optional pairing, and off-target detection.
"""

import ap
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import re
import argparse
import os
import sys

excavate_banner = r"""

    _______  ___________ _    _____  ____________     __  ________
   / ____/ |/ / ____/   | |  / /   |/_  __/ ____/    / / / /_  __/
  / __/  |   / /   / /| | | / / /| | / / / __/______/ /_/ / / /   
 / /___ /   / /___/ ___ | |/ / ___ |/ / / /__/_____/ __  / / /    
/_____//_/|_\____/_/  |_|___/_/  |_/_/ /_____/    /_/ /_/ /_/     
    
    ⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⣀⣀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⢀⣀⡿⠿⠿⠿⠿⠿⠿⢿⣀⣀⣀⣀⣀⡀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠸⠿⣇⣀⣀⣀⣀⣀⣀⣸⠿⢿⣿⣿⣿⡇⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠻⠿⠿⠿⠿⠿⣿⣿⣀⡸⠿⢿⣿⡇⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣿⣿⣿⣧⣤⡼⠿⢧⣤⡀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣿⣿⣿⣿⠛⢻⣿⡇⠀⢸⣿⡇
    ⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣿⣿⣿⣿⠛⠛⠀⢸⣿⡇⠀⢸⣿⡇
    ⠀⠀⠀⠀⠀⠀⢠⣤⣿⣿⣿⣿⠛⠛⠀⠀⠀⢸⣿⡇⠀⢸⣿⡇
    ⠀⠀⠀⠀⢰⣶⣾⣿⣿⣿⠛ ⠀⠀⠀⠀⠀⠈⠛⢳⣶⡞⠛⠁⠀
    ⢰⣶⡎⠉⢹⣿⡏⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⢸⣿⣷⣶⡎⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                                                     
    """
def add_common_args(parser):

    parser.add_argument(
        '--pairing-method',
        nargs='?',
        const='r',
        choices=['r', 'fp', 't'],
        type=str,
        help="Enable pairing of gRNA to output a dual-guide library. Default is 'r': pairing all guides that target different SNPs together. 'fp': pairing guides about a fixed point. 't': tiled pairing"
    )
    
    parser.add_argument(
        "-f",
        "--fixed-points-list",
        type=str,
        help="1 or more fixed points, comma-separated (genomic coordinate without chr#, eg: 11989251,12002042,...]) in your locus to pair guides around."
    )
        
    """Add arguments common to multiple subcommands"""
    parser.add_argument(
       "-o", "--output_dir",
        required=True,
        type=str,
        help="Output directory. Folder name if you are already in the excavate folder.",
    )
        

        
def add_generate_parser(subparsers):
    
    """Add the generate subcommand parser"""
    gen_parser = subparsers.add_parser(
        "generate", 
        help="Generate complete gRNA library using inputted variant data"
    )

    gen_parser.add_argument("vcf_file", type=str, help="path to the vcf file, comma-separated if more than one (cell_line.vcf.gz,1000genomes.vcf.gz). If using a phased cell-line vcf, put that first.")
    
    gen_parser.add_argument("var_type", type=str, help="'cell-line' or 'population', with a number and comma-separated if more than one (eg: cell-line,population1,population2)")
    
    gen_parser.add_argument("fa_file", type=str, help="path to the chromosome fasta file for your locus of interest (eg: ch1sequence.fasta)")
    
    gen_parser.add_argument("genome_fa", type=str, help="path to the whole genome fasta file for your organism")
    
    gen_parser.add_argument("locus", type=str, help="chr#:start-end")
    
    gen_parser.add_argument(
        "--af-threshold", 
        type=str, 
        nargs="?",
        default='0.1',
        const='0.1', 
        help="allele frequency threshold between 0 and 1 (default = 0.1)"
    )

    gen_parser.add_argument( #this needs to be an option, not a flag
        "--cas",
        type=str,
        help="One of SpCas9, SpCas9_NG, enAsCas12a, or SaCas9",
    )
    
    gen_parser.add_argument( #this needs to be an option, not a flag
        "--pam-list",
        type=str,
        help="PAM sequences, if not using one of the supported Cas-species. Comma-separated. Use IUPAC codes.",
    )
    
    gen_parser.add_argument( #this needs to be an option, not a flag
        "--orient",
        choices=['3prime', '5prime'],
        type=str,
        help="PAM orientation",
    )
    
    gen_parser.add_argument(
        "-g", 
        "--guide-length",
        type=int, 
        nargs="?",
        default=20,
        const=20, 
        help="guide length in base pairs (defaut = 20)"
    )
    
    gen_parser.add_argument(
        "-m",
        "--max_snppos_in_protospacer", 
        type=int, 
        nargs="?",
        default=10,
        const=10, 
        help="maximum distance in bp of SNP from PAM sequence (default = 10)"
    )
                        
    gen_parser.add_argument(
        "--off-targets",
        action='store_true',
        help="Include for off-targets analysis (counting exact matches in genome, and 1 bp mismatches in chromosome of interest).",
    )
                        
    gen_parser.add_argument(
        "--split-phased",
        action='store_true',
        help="Enable splitting of gRNA libraries by cell-line phasing.",
    )

    gen_parser.add_argument(
        "--summary",
        action='store_true',
        help="Enable to output a summary table for each gRNA library.",
    )

    gen_parser.add_argument(
        "--per-vcf",
        action='store_true',
        help="Enable to save single-gRNA libraries for each VCF file, split by allele.",
    )

    add_common_args(gen_parser)
    return gen_parser


def add_pair_parser(subparsers):
    """Add the pair subcommand parser"""
    pair_parser = subparsers.add_parser(
        "pair",
        help="Pair your single-gRNA library"
    )
    pair_parser.add_argument(
        "-i", "--input-library", 
        required=True,
        type=str,
        help="Path to input single-gRNA library file"
    )

    add_common_args(pair_parser)
    return pair_parser



def initialize_cas_obj(args):
    
    # Initialize cas_obj based on inputted Cas enzyme or custom PAM and orientation
    if args.cas is not None and args.pam_list is not None:
        raise ValueError('Both Cas and PAM-list given, please only input one of these.')
        
    elif args.cas is not None:
        cas_name = args.cas
        if hasattr(ap, cas_name):
            cas_obj = getattr(ap, cas_name)
        else:
            raise ValueError(f"Invalid Cas enzyme name: {cas_name}, choose one of SpCas9, SpCas9-NG, enAsCas12a, or SaCas9. Otherwise, input custom PAMs")

    elif args.pam_list is not None:
        if args.orient is not None:
            pam_list = args.pam_list.split(',')
            cas_obj = ap.create_custom_cas_obj(pam_list, args.orient)
        else:
            raise ValueError("Custom PAM orientation not provided. Input --orient 5prime or 3prime")
            
    else:
        raise ValueError("Neither Cas nor custom PAM provided. Please input one of these.")

    return cas_obj

def load_variant_data(vcf_file, var_type, locus, af_threshold):

    # check if vcf file path exists
    for f in vcf_file:
        if not os.path.exists(f):
            raise FileNotFoundError(f"VCF file '{f}' does not exist.")

    # check if var_types are valid
    valid_types = {'cell-line', 'population'}

    def get_chars_before_number(text): # function to store only the letter chars from var_type, to aid checking if the input is valid
        for i, char in enumerate(text):
            if char.isdigit():
                return text[:i]
        return text
        
    for vt in var_type:
        vt_letters = get_chars_before_number(vt)
        if vt_letters not in valid_types:
            raise ValueError(f"Invalid variant type: {vt}. Must be one of {valid_types}.")

    # check if all inputs are the same length
    if not (len(vcf_file) == len(var_type)):
        raise ValueError("Arguments vcf_file and var_type must have the same number of comma-separated values.")
        
    # check if af threshold is a valid float
    try:
        af_threshold = float(af_threshold)
        if not (0 <= af_threshold <= 1):
            raise ValueError("Allele frequency threshold must be between 0 and 1.")
    except ValueError:
        raise ValueError("Invalid allele frequency threshold. Must be a number between 0 and 1.")
            
    # executing create_gens to create dataframe of SNP genotypes from VCF file data

    gens_dict = {}  # Dictionary to store gens dataframes
    for i in range(len(vcf_file)): # For each vcf file provided
        key = f"gens_df_{var_type[i]}"  # Create a unique key for each iteration
        gens_dict[key] = ap.create_gens(vcf_file[i], locus, var_type[i], af_threshold)
        
    for i in range(len(vcf_file)):
        key = f"gens_df_{var_type[i]}" # Create a unique key for each iteration
        try:
            gens_dict[key] = ap.create_gens(vcf_file[i], locus, var_type[i], af_threshold)
        except Exception as e:
            raise RuntimeError(f"Failed to create gens data for file {vcf_file[i]} with variant type {var_type[i]}: {e}")

    return gens_dict

def run_off_targets(all_guides_unique, genome_fa, cas_obj, chseq):

    def prompt_continue():
        while True:
            response = input("This step may take several hours. Do you want to continue? [y/n]: ").strip().lower()
            if response in ('y', 'yes'):
                return True
            elif response in ('n', 'no'):
                return False
            else:
                print("Please respond with 'y' or 'n'.")        
                
    # Prompt user
    if prompt_continue():
        
        print('Off-targets analysis has begun, but may take hours to finish. To allow it to run without terminating, prevent your computer from sleeping. For Mac, go to Display > Advanced > Prevent automatic sleeping... Additionally, please make sure your computer does not shut down or run out of power')

        try:
            all_guides_exact_counts = ap.count_exact_matches(all_guides_unique, genome_fasta_path=genome_fa, cas_parameters=cas_obj)
            all_gudies_exact_1bp_counts = ap.one_mismatch(all_guides_exact_counts, chseq)
            all_guides_final = all_gudies_exact_1bp_counts
            print('Off-target analysis completed!')
        except Exception as e:
            raise RuntimeError(f"Failed to peform off-targets analysis: {e}")
        
    else:
        all_guides_final = all_guides_unique
        print('Okay! Exiting out of off-target analysis...')
        
    return all_guides_final
        
def apply_pairing(all_guides, method, fixed_points=None):
    if method == 'r':
        return ap.random_pair(all_guides)
    elif method == 'fp':
        return ap.fixed_point_pair(all_guides, fixed_points)
    elif method == 't':
        return ap.tiled_pair(all_guides)
    else:
        raise ValueError(f"Unsupported pairing method: {method}")

def run_generate(args):

    print(excavate_banner)
    
    # get output directory
    outdir = args.output_dir

    cas_obj = initialize_cas_obj(args)
    print(f'Initialized Cas object: {cas_obj.name}')

    # Use inputted VCF files to generate dataframe of SNPs to find gRNA for.

    # Extract all args from args
    
    vcf_file = args.vcf_file.split(',')
    var_type = args.var_type.split(',')
    
    #checking validity of locus arg, and saving start and end for downstream uses (mainly verifying if inputted fixed-points are within the locus)
    locus = args.locus
    region_match = re.match(r"chr\w+:(\d+)-(\d+)", args.locus)
    if not region_match:
        raise ValueError("Locus format must be chr:start-end")
    start = int(region_match.group(1))
    end = int(region_match.group(2))

    af_threshold = args.af_threshold
    
    # create list of snps according to criteria, for each given vcf file and var type
    gens_dict = load_variant_data(vcf_file, var_type, locus, af_threshold)
    print(f'Loaded variant data from {vcf_file}')

    # Make chr Seq object using chr fasta - this is not in a function, shoud it be? maybe
    
    fa_file = args.fa_file
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"Chromosome fasta file '{fa_file}' does not exist.")
    else:
        try:
            chseq = ap.makeseq(fa_file)
        except Exception as e:
            raise RuntimeError(f"Failed to load chromosome sequence from file {fa_file}: {e}")
    print(f'Loaded chromosome sequence from {fa_file}')

    # Create allelic sequences for each VCF file, with each of the present versions of the SNPs
    
    seq_dict = {}
    try:
        for i in gens_dict: # for each set of snps
            key = f"{i}_allele1_seq"
            seq_dict[key] = ap.getaltseq(gens_dict[i], chseq, snpform = 'allele1')
            key = f"{i}_allele2_seq"
            seq_dict[key] = ap.getaltseq(gens_dict[i], chseq, snpform = 'allele2')
    except Exception as e:
        raise RuntimeError(f"Failed to create allelic sequences from {gens_dict[i]}: {e}")
    
    # Find guides using gens_dict[i remove the last 8 chars, to give name of corresponding gens df], for each seq.

    # check if args for find_guides are valid
    max_snppos_in_protospacer = args.max_snppos_in_protospacer
    guide_length = args.guide_length
    
    if not (1 <= guide_length):
        raise ValueError("guide length must be at least 1.")
    
    if not (1 <= max_snppos_in_protospacer <= guide_length):
        raise ValueError("maximum SNP pos in protospacer must be between 1 and guide length.")

    # run find guides with valid args
    guides_dict = {}
    for i in seq_dict:
        key = f"{i}_guides"
        guides_dict[key] = ap.find_guides(gens_dict[str(i)[:-12]], seq_dict[i], cas_obj, max_snppos_in_protospacer, guide_length)

    # save individual guides dfs in separate folder in output directory if --per-vcf included
    if args.per_vcf:
        
        os.makedirs(os.path.join(outdir, "single_gRNA_libraries_per_vcf"), exist_ok=True)
        for lib in guides_dict:
            single_gRNA_libs_per_vcf_output = os.path.join(outdir, "single_gRNA_libraries_per_vcf", f"{lib}.csv")
            guides_dict[lib].to_csv(single_gRNA_libs_per_vcf_output, index=False)
    
        print(f'Saved single-gRNA libraries for each VCF file provided, split by allele, in {outdir}/single_gRNA_libraries_per_vcf.')

    # combine all gRNA libraries and add annotations about variant presence and alt allele frequencies
    print('Combining all single-gRNA libraries...')
    all_guides = pd.concat(list(guides_dict.values()), ignore_index=True)
    all_guides_unique = ap.all_guides_var_info(gens_dict, all_guides)
    all_guides_unique = all_guides_unique.reset_index(drop = True)

    # Off-target analysis, if enabled
    
    if args.off_targets:
        print('Preparing to run off-targets analysis...')
        genome_fa = args.genome_fa
        if not os.path.exists(genome_fa):
            raise FileNotFoundError(f"Whole genome fasta file '{genome_fa}' does not exist.")
        else:
            all_guides_final = run_off_targets(all_guides_unique, genome_fa, cas_obj, chseq)
        # run_off_targets() which returns all_guides_final
    else:
        all_guides_final = all_guides_unique

    # Save final all_guides_unique df with or without off-targets analysis in outdir
    
    # make sure outdir exists (although we already did this above, can maybe remove)
    os.makedirs(outdir, exist_ok=True)

    # save all_guides df in output directory as a csv file
    all_guides_final_output = os.path.join(outdir, "all_guides.csv")
    all_guides_final.to_csv(all_guides_final_output, index=False)

    print(f'Combined and annotated single-gRNA library saved in {outdir}')

    # Split by phasing, if enabled

    # check validity of fixed-points-lists if enabled. Checking here because pairing happens in each if block sepately. Do not want to check validity multiple times.
    # first, check if fixed-points given if pair == fp
    fixed_points_list = args.fixed_points_list
    
    if args.pairing_method == 'fp':
        if fixed_points_list is None:
            raise ValueError("Fixed point pairing enabled, please input 1 or more fixed points using --fixed-points-list")
        else:
            fixed_points_list = fixed_points_list.split(',')

    # second, check if points are valid (int, and fall between locus start and end)
        for point in fixed_points_list:
            try:
                point_int = int(point)
            except ValueError:
                raise ValueError(f"Invalid fixed point: {point} is not an integer")
            if not (start < point_int < end):
                raise ValueError(f"Fixed point {point} is out of locus range {start}-{end}")


    # if all is gucci, go ahead with splitting by phasing and pairing

    if args.split_phased:
    # split
        phased_vcf = vcf_file[0]
        try:
            split_list = ap.split_phased(all_guides_final, phased_vcf, locus)
            all_guides_allele1 = split_list[0]
            all_guides_allele2 = split_list[1]
            
            all_guides_allele1_output = os.path.join(outdir, "all_guides_allele1.csv")
            all_guides_allele2_output = os.path.join(outdir, "all_guides_allele2.csv")

            all_guides_allele1.to_csv(all_guides_allele1_output, index=False)
            all_guides_allele2.to_csv(all_guides_allele2_output, index=False)

            print(f'Single-gRNA library split by phasing, saved in {outdir}')
            
        except Exception as e:
            raise RuntimeError(f"Failed to split all_guides dataframe by phased alleles: {e}")
            
        # if pairing enabled, pair each allele's guides according to specified method    
        if args.pairing_method:
            all_guides_allele1_paired = apply_pairing(all_guides_allele1, args.pairing_method, fixed_points_list)
            all_guides_allele2_paired = apply_pairing(all_guides_allele2, args.pairing_method, fixed_points_list)

            all_guides_allele1_paired_output = os.path.join(outdir, "all_guides_allele1_paired.csv")
            all_guides_allele2_paired_output = os.path.join(outdir, "all_guides_allele2_paired.csv")

            all_guides_allele1_paired.to_csv(all_guides_allele1_paired_output, index=False)
            all_guides_allele2_paired.to_csv(all_guides_allele2_paired_output, index=False)

            print(f'Paired dual-guide libraries for each allele saved in {outdir}')
            
    # if split-phased disabled, don't split, perform pairing on all_guides_final        
    else:
        if args.pairing_method:
            all_guides_paired = apply_pairing(all_guides_final, args.pairing_method, fixed_points_list)
            
            all_guides_paired_output = os.path.join(outdir, "all_guides_paired.csv")
            all_guides_paired.to_csv(all_guides_paired_output, index=False)

            print(f'Paired dual-guide library saved in {outdir}')
            
    # Generate summary table if --summary enabled

    if args.summary:
        if args.split_phased:
            all_guides_allele1_summary = ap.targetable_vars(all_guides_allele1)
            all_guides_allele2_summary = ap.targetable_vars(all_guides_allele2)

            all_guides_allele1_summary_output = os.path.join(outdir, "all_guides_allele1_summary.csv")
            all_guides_allele2_summary_output = os.path.join(outdir, "all_guides_allele2_summary.csv")

            all_guides_allele1_summary.to_csv(all_guides_allele1_summary_output, index=False)
            all_guides_allele2_summary.to_csv(all_guides_allele2_summary_output, index=False)

            print('Summary tables saved')
            
        else:
            all_guides_summary = ap.targetable_vars(all_guides_final)
            
            all_guides_summary_output = os.path.join(outdir, "all_guides_summary.csv")
            all_guides_summary.to_csv(all_guides_summary_output, index=False)
            
            print('Summary table saved')


    print('EXCAVATE has finished running! Enjoy your gRNA libraries! :)')
    
    return

def run_pairing(args):

    print(excavate_banner)
    
    # get output directory
    outdir = args.output_dir
    
    #sort input library by values of SNP position and guide start
    input_library_path = args.input_library

    if not os.path.exists(input_library_path):
        raise FileNotFoundError(f"input library file '{input_library_path}' does not exist.")

    all_guides = pd.read_csv(input_library_path)
    all_guides = all_guides.sort_values(by=["SNP position","start"]).reset_index(drop=True)

    #check if fp
    fixed_points_list = args.fixed_points_list
    
    if args.pairing_method == 'fp':
        if fixed_points_list is None:
            raise ValueError("Fixed point pairing enabled, please input 1 or more fixed points using --fixed-points-list")
        else:
            fixed_points_list = fixed_points_list.split(',')

    # second, check if points are valid (int, and fall between locus start and end)
        for point in fixed_points_list:
            try:
                point_int = int(point)
            except ValueError:
                raise ValueError(f"Invalid fixed point: {point} is not an integer")
            
    # if pairing enabled, pair each allele's guides according to specified method    
    all_guides_paired = apply_pairing(all_guides, args.pairing_method, fixed_points_list)
    if all_guides_paired.empty:
        raise ValueError("No pairs made. Check if fixed-points fall within library's genomic region")

    os.makedirs(outdir, exist_ok=True)
    all_guides_paired_output = os.path.join(outdir, "all_guides_paired.csv")
    
    all_guides_paired.to_csv(all_guides_paired_output, index=False)

    print(f'Paired dual-guide library saved in {outdir}')

    return
        

def main():
    
    """Main entry point for the gRNA pipeline"""
    parser = argparse.ArgumentParser(
        description="ExCAVaTE: a pipeline to identify targetable genomic variants and generate allele-specific single and paired-gRNA libraries",
        prog="excavate"
    )
    
    # Add version argument
    parser.add_argument(
        "--version", 
        action="version", 
        version="%(prog)s 1.0.0"
    )
    
    # Create subparsers
    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="Available pipeline commands",
        help="Use 'excavate <command> --help' for command-specific help"
    )
    
    # Add subcommand parsers
    gen_parser = add_generate_parser(subparsers)
    pair_parser = add_pair_parser(subparsers)

     # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()

    # Handle case where subcommand is provided but no arguments
    if not hasattr(args, 'output_dir'):
        if args.command == "generate":
            gen_parser.print_help()
        elif args.command == "pair":
            pair_parser.print_help()
        else:
            parser.print_help()
        sys.exit(0)

    # Execute appropriate command
    if args.command == "generate":
        run_generate(args)
    elif args.command == "pair":
        run_pairing(args)

if __name__ == "__main__":
    main()