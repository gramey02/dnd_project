import numpy as np
import pandas as pd
import re
import regex
import subprocess
from io import StringIO
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
import multiprocessing as mp
from pyfaidx import Fasta

"""
Defining the Cas class. And defininh all functions to create Cas objects, and creating some set Cas objects 
"""
class Cas:
    def __init__(self, name, pam_three_prime, pam_five_prime, pam_length,
                 exclusion_positions, is_five_prime, is_three_prime, is_multi_pam):
        self.name = name
        self.pam_three_prime = pam_three_prime
        self.pam_five_prime = pam_five_prime
        self.pam_length = pam_length
        self.exclusion_positions = exclusion_positions
        self.is_five_prime = is_five_prime
        self.is_three_prime = is_three_prime
        self.is_multi_pam = is_multi_pam

SpCas9 = Cas(
    name='SpCas9',
    pam_three_prime=r'(?P<p0>[ACTG]GG)',
    pam_five_prime=r'(?P<p0>CC[ACTG])',
    pam_length=3,
    exclusion_positions=[[-1]],
    is_five_prime=False,
    is_three_prime=True,
    is_multi_pam=False
)

enAsCas12a = Cas(
    name='enAsCas12a',
    pam_three_prime=r'(?P<p0>[ACGT][GA]AA)|(?P<p1>[TGC]AA[TGC])|(?P<p2>[TGC]A[TC]A)|(?P<p3>GG[ACGT][TG])',
    pam_five_prime=r'(?P<p0>TT[CT][ACGT])|(?P<p1>[ACG]TT[ACG])|(?P<p2>T[AG]T[ACG])|(?P<p3>[AC][ACGT]CC)',
    pam_length=4,
    exclusion_positions=[[-1, 2], [-1, -4], [-1], [-3, 2]],
    is_five_prime=True,
    is_three_prime=False,
    is_multi_pam=True
)

SpCas9_NG = Cas(
    name='SpCas9-NG',
    pam_three_prime= r'(?P<p0>[ACTG]G)',
    pam_five_prime=r'(?P<p0>C[ACTG])',
    pam_length=2,
    exclusion_positions=[[-1]],
    is_five_prime=False,
    is_three_prime=True,
    is_multi_pam=False
)

SaCas9 = Cas(
    name='SaCas9',
    pam_three_prime=r'(?P<p0>[ACTG][ACTG]G[AG][AG]T)',
    pam_five_prime=r'(?P<p0>A[TC][TC]C[ACTG][ACTG])',
    pam_length=6,
    exclusion_positions=[[-1, -2]],
    is_five_prime=False,
    is_three_prime=True,
    is_multi_pam=False
)

"""
function to output a regex pattern from IUPAC PAM code.
"""
def create_regex_pam(custom_pam_list):
    
    replace_dict = {
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B':'[CGT]',
    'D': '[AGT]',
    'H':'[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]'
    }

    custom_pam_list_regex = []
    for pam in custom_pam_list:
        for char in pam:
            if char in replace_dict.keys():
                pam = pam.replace(char, replace_dict[char])
            else:
                pass
        custom_pam_list_regex.append(pam)
                    
    return custom_pam_list_regex

"""
function to create the exclusion position set for a given custom PAM list, and orientation
"""
def create_exclusion_pos_set(custom_pam_list, orient):
    
    exclusion_codes = 'BDHVN'
    all_indices = []
    pam_len = len(custom_pam_list[0])

    for pam in custom_pam_list:
        pamwise_indices = []
        for index, code in enumerate(pam):
            if code in exclusion_codes:
                pamwise_indices.append(index)
        all_indices.append(pamwise_indices)
    
    exclusion_pos = []
    if orient == '5prime':
        for pamwise_indices in all_indices:
            pamwise_exclusion_pos = []
            for index in pamwise_indices:
                pos = index - pam_len
                pamwise_exclusion_pos.append(pos)
            exclusion_pos.append(pamwise_exclusion_pos)
    else:
        for pamwise_indices in all_indices:
            pamwise_exclusion_pos = []            
            for index in pamwise_indices:
                pos = - 1 - index
                pamwise_exclusion_pos.append(pos)
            exclusion_pos.append(pamwise_exclusion_pos)
            
    return exclusion_pos

"""
function to create a custom cas_obj
"""
def create_custom_cas_obj(custom_pam, orient):
    #custom_pam can have multiple values (multi-PAM)

    pam_len = len(custom_pam[0])

    custom_pam_rc = []
    for pam in custom_pam:
        pam_seq = Seq(pam)
        pam_seq_rc = pam_seq.reverse_complement()
        custom_pam_rc.append(str(pam_seq_rc))
    
    custom_pam_regex = '|'.join(f"(?P<p{i}>{pam})" for i, pam in enumerate(create_regex_pam(custom_pam)))
    custom_pam_rc_regex = '|'.join(f"(?P<p{i}>{pam})" for i, pam in enumerate(create_regex_pam(custom_pam_rc)))
    
    if orient == '5prime':
        pam_five_prime = custom_pam_regex
        pam_three_prime = custom_pam_rc_regex
        is_five_prime=True
        is_three_prime=False
    else:
        pam_five_prime = custom_pam_rc_regex
        pam_three_prime = custom_pam_regex
        is_five_prime=False
        is_three_prime=True

    if len(custom_pam) > 1:
        is_multi_pam = True
    else:
        is_multi_pam = False

    exclusion_positions = create_exclusion_pos_set(custom_pam, orient)
    
    custom_cas = Cas(
        name='custom_cas',
        pam_length=pam_len,
        pam_three_prime=pam_three_prime,
        pam_five_prime=pam_five_prime,
        exclusion_positions=exclusion_positions,
        is_five_prime=is_five_prime,
        is_three_prime=is_three_prime,
        is_multi_pam=is_multi_pam
    )

    return custom_cas


# ESSNTIAL FUNCTIONS: PREPROCESSING VCF FILES AND FASTA FILES INPUTTED.

"""
Import fasta file for required chromosome
"""
def makeseq(fa_filename):
    file_path = fa_filename
    with open(file_path, 'r') as fasta_file:
        # Parse the FASTA file
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence_id = record.id
            sequence = record.seq
    chseq = record.seq
    return chseq


"""
Returns a table of genotypes at each SNP location, with SNPs filtered by genomic region, heterozygosity in a cell line, or allele frequency in a population database
Arguments:
    vcf_filename: string filename
    genomic_region: chr#:startpos-endpos
    variants_db: 'population' or 'cell-line' to handle filtering of SNPs differently
    af_threshold: default 0.3, select SNPs with minor allele frequency above this threshold 
"""
def create_gens(vcf_file, locus, variants_db, af_threshold = 0.1): #for snps, not indels rn

    #get chrom formatting in vcf file
    cmd = f"bcftools view {vcf_file} | grep -v '^#' | head -n 1 | cut -f1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    chrom_formatting = result.stdout.strip()

    if chrom_formatting.startswith('chr'):
        pass
    else:
        locus = locus[3:]
        
    # handle population and cell-line variant databases differently because population will use allele frequency
    if variants_db.startswith('population'):
        bcl_view = subprocess.Popen(f'bcftools view -v snps -g ^miss -q {af_threshold}:minor -r {locus} {vcf_file} -Ou | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n"', shell=True, stdout=subprocess.PIPE)
        bcl_view.wait()
        col_names_sim = ["chrom", "pos", "snp_id", "ref", "alt", "alt AF"]

    elif variants_db.startswith('cell-line'):
        bcl_view = subprocess.Popen(f'bcftools view -v snps -g ^miss -g het -r {locus} {vcf_file} -Ou | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"', shell=True, stdout=subprocess.PIPE)
        bcl_view.wait()
        col_names_sim = ["chrom", "pos", "snp_id", "ref", "alt", 'genotype']

    else:
        raise ValueError("var_types must be either 'cell-line' or 'population' optionally followed by numbers if inputting mutiple of each type")
    
    gens = pd.read_csv(
            StringIO(bcl_view.communicate()[0].decode("utf-8")),
            sep="\t",
            header=None)
    gens.columns = col_names_sim

    #add another column to gens to store the present alleles acccording to the genotype. Hence, for cell-line vcfs, if the file is phased, this is information about which alleles are phased together. If file is not phased, this helps store the alleles present in the cell-line in case any of the snps are multi-allelic, and the ref version isn't present. This prevents gRNA being made with a version of the snp that isn't present in the cell-line.
    present_alleles_list = []
    for i in range(len(gens)):
        alleles = gens['ref'][i] + '(ref)' + ',' + gens['alt'][i]
        alleles = alleles.split(',')
        for a in range(len(alleles)):
            if a > 0:
                alleles[a] = alleles[a] + '(alt)'
            else:
                continue

        if variants_db.startswith('cell-line'):
            gt = gens['genotype'][i]
            gt = regex.split(r'[\/|]', gt)
            allele_index1 = int(gt[0])
            allele_index2 = int(gt[1])
        else:
            #if its a population gensdf, I am choosing to not deal with multi-allelic snps at the moment. simpy use ref and the first alt allele as the present alleles.
            #want to have a uniform present alleles column in cell-line and population gensdfs so we can extract allelic information downstream in the same way for everything
            allele_index1 = 0
            allele_index2 = 1

        #using set allele_indices, extract the corresponding alleles to generate the present_alleles data
        present_alleles = (alleles[allele_index1], alleles[allele_index2])
        present_alleles_list.append(present_alleles)
        
    gens['present alleles'] = present_alleles_list

    return gens


"""
Returns an immutable sequence with either ref or alt forms of SNP targets of interest
Arguments:
    vcf : string filename
    refseq : reference sequence to be mutated
    snpform = 'ref' or 'alt'
    dbsnp = True or False (default False) to handle case of using dbsnp variant file input
"""
def getaltseq(gens_df, refseq, snpform = 'allele1'):
    mutable_seq = MutableSeq(refseq)
    gens_df['pos'] = gens_df['pos'].astype(int)
    for i in range(len(gens_df)):
        if snpform == 'allele1':
            #picks out ith eentry in present_alleles col (eg: (T(ref), C(alt))), 0th entry of the tuple (T(ref)), and 0th element of that string (T) to get the present allele in allele1
            present_allele = gens_df['present alleles'][i][0][0]
            mutable_seq[gens_df['pos'][i]-1] = present_allele
        elif snpform == 'allele2':
            #same thing for allele2, but 1th entry in the tuple to get the allele in allele2
            present_allele = gens_df['present alleles'][i][1][0]
            mutable_seq[gens_df['pos'][i]-1] = present_allele
            
    immutable_seq = Seq(mutable_seq)
    return immutable_seq


# MEAT OF THE PIPELINE: FIND GUIDES FUNCTION

"""
Find guides in given sequence around given loci for particular Cas species
"""
def find_guides(snplist, sequence, cas_obj, max_snp_pos_in_protospacer, guide_len):

    pam_len = cas_obj.pam_length
    
    guidelist = []
    PAMlist = []
    offtarget_counts_list = []
    chrom_list = []
    start_list = []
    end_list = []
    SNPcoordinates = []
    rsIDs = []
    variation = []
    strand = []
    SNPallele = []
    guidename = []
    pos_in_protospacer_or_pam = [] #1,2,3... for pos in protospacer, or -1,-2,-3 for pos in PAM starting from N
        
    for i in range(len(snplist)):
    
        snppos_index = int(snplist['pos'][i]-1)
        snpregion = str(sequence[snppos_index-(guide_len+(pam_len-1)):snppos_index+(guide_len+pam_len)]) #snp is at sequence[pos-1]
        gc_index_conversion_num = snppos_index-(guide_len+(pam_len-1)) + 1 # index + gc_index_conversion_num = genomic coordinate index in string, +1 = genomic coordinate
    
        #find downstream pams on plus strand
        plus_starts = []
        plus_matches = []
        plus_matched_group_list = []
        
        for pam in regex.finditer(cas_obj.pam_three_prime, snpregion, regex.IGNORECASE, overlapped=True):        
            plus_starts.append(pam.start())
            plus_matches.append(pam.group(0))
            matched_group = None
            for name in pam.groupdict():
                if pam.group(name):
                    matched_group = name
                    break
            plus_matched_group_list.append(matched_group)
            
        #find upstream pams on minus strand
        minus_starts = []
        minus_matches = []
        minus_matched_group_list = []
        
        for pam in regex.finditer(cas_obj.pam_five_prime, snpregion, regex.IGNORECASE, overlapped=True):
            minus_starts.append(pam.start())
            minus_matches.append(pam.group(0))
            matched_group = None
            for name in pam.groupdict():
                if pam.group(name):
                    matched_group = name
                    break
            minus_matched_group_list.append(matched_group)
       
        snppos_in_snpregion = guide_len+pam_len-1
    
        chrom_num = str(re.findall(r'\d+', str(snplist['chrom'][0]))[0])
    
        #find guides on plus strand
        for j in range(len(plus_starts)):
            start = plus_starts[j]
            PAM = plus_matches[j]
            matched_group = plus_matched_group_list[j]
            if ((start >= guide_len) and (start <= snppos_in_snpregion + max_snp_pos_in_protospacer)): #all the way up to the SNP being in the last position on the PAM - because that is always concrete, to N of PAM being max pos away from snp.
    
                #compute position of SNP in protospacer or PAM, and discard guide if SNP is in a position that should be excluded
                pos_in_protospacer_or_pam_num = start-snppos_in_snpregion
                if pos_in_protospacer_or_pam_num <=0: #if position is in PAM, subtract 1 to follow convention of -1,-2,-3 being the positions of the PAM
                    pos_in_protospacer_or_pam_num = pos_in_protospacer_or_pam_num-1
                else:
                    pass
    
                exclusion_set = cas_obj.exclusion_positions[int(matched_group[-1])]
                if pos_in_protospacer_or_pam_num in exclusion_set:
                    flag = True
                else:
                    flag = False
                if flag:
                    continue
                else:
                    pos_in_protospacer_or_pam.append(pos_in_protospacer_or_pam_num)
    
                #saving PAM and guide, differently handled based on 5 prime or 3 prime PAM
                if cas_obj.is_five_prime:
                    PAMlist.append(str(Seq(PAM).reverse_complement()))
                    guide = Seq(snpregion[start-guide_len: start])
                    guidelist.append(str(guide.reverse_complement()))
                else:
                    PAMlist.append(str(PAM))
                    guide = snpregion[start-guide_len: start]
                    guidelist.append(str(guide))
    
                chrom = ('chr' + chrom_num)
                chrom_list.append(chrom)
                
                guide_start = start - guide_len + gc_index_conversion_num #pam_start - guide_len + conversion_num
                start_list.append(guide_start)
    
                guide_end = guide_start + guide_len - 1 # end of the guide, inclusive
                end_list.append(guide_end)
                
                SNPcoordinates.append(snplist['pos'][i])
    
                if 'snp_id' in snplist.columns:
                    rsIDs.append(snplist['snp_id'][i])
                else:
                    rsIDs.append('None')
    
                variation.append((snplist['ref'][i] + '>' + snplist['alt'][i]))

                if cas_obj.is_five_prime: #save guides as minus guides if cas PAM is five-prime, because that's how the output makes sense
                    guide_strand = "-"
                else:
                    guide_strand = "+"
                
                strand.append(guide_strand)
                    
                presentallele = sequence[snppos_index]
                if presentallele == snplist['ref'][i]:
                    SNPallele.append(presentallele + ' (ref)') #go to the specific snp position in the given sequence and print that
                else:
                    SNPallele.append(presentallele + ' (alt)') #go to the specific snp position in the given sequence and print that
    
                guidename.append('chr' + chrom_num + '_' + str(guide_start) + '_' + guide_strand + '_' + presentallele + '_' + cas_obj.name + '_' + str(guide_len) + 'nt')

            else:
                pass
    
        #find guides on minus strand
        for j in range(len(minus_starts)):
            start = minus_starts[j]
            PAM = minus_matches[j]
            matched_group = minus_matched_group_list[j]
    
            if (start >= (snppos_in_snpregion - max_snp_pos_in_protospacer) and (start <= snppos_in_snpregion)): #all the way up to the SNP being in the last position on the PAM - because that is always concrete, to N of PAM being max pos away from snp
    
                pos_in_protospacer_or_pam_num = -1*(start + pam_len-1 - snppos_in_snpregion) #such that 'N' is -1, protospacer positions are positive integers, and PAM positions are negative
                if pos_in_protospacer_or_pam_num <=0:
                    pos_in_protospacer_or_pam_num = pos_in_protospacer_or_pam_num-1
                else:
                    pass
    
                exclusion_set = cas_obj.exclusion_positions[int(matched_group[-1])]
                if pos_in_protospacer_or_pam_num in exclusion_set:
                    flag = True
                else:
                    flag = False
                if flag:
                    continue
                else:
                    pos_in_protospacer_or_pam.append(pos_in_protospacer_or_pam_num)
                
                if cas_obj.is_five_prime:
                    PAMlist.append(str(PAM))
                    guide = snpregion[start+pam_len : start+pam_len+guide_len]
                    guidelist.append(str(guide))
                    
                else:
                    PAMlist.append(str(Seq(PAM).reverse_complement()))
                    guide = Seq(snpregion[start+pam_len : start+pam_len+guide_len])
                    guidelist.append(str(guide.reverse_complement()))
    
                chrom = ('chr' + chrom_num)
                chrom_list.append(chrom)
                
                guide_start = start + pam_len + gc_index_conversion_num
                start_list.append(guide_start)
    
                guide_end = guide_start + guide_len - 1 # end of the guide, inclusive
                end_list.append(guide_end)
                
                # guide_coords = (guide_start, (guide_start + guide_len-1)) #start inclusive
                # guide_coords_list.append(guide_coords)
                
                SNPcoordinates.append(snplist['pos'][i])
    
                if 'snp_id' in snplist.columns:
                    rsIDs.append(snplist['snp_id'][i])
                else:
                    rsIDs.append('None')
    
                variation.append((snplist['ref'][i] + '>' + snplist['alt'][i]))
    
                if cas_obj.is_five_prime: #save guides as minus guides if cas PAM is five-prime, because that's how the output makes sense
                    guide_strand = '+'
                else:
                    guide_strand = '-'
                    
                strand.append(guide_strand)
                    
                presentallele = sequence[snppos_index]
                if presentallele == snplist['ref'][i]:
                    SNPallele.append(presentallele + ' (ref)') #go to the specific snp position in the given sequence and print that
                else:
                    SNPallele.append(presentallele + ' (alt)') #go to the specific snp position in the given sequence and print that
    
                guidename.append('chr' + chrom_num + '_' + str(guide_start) + '_' + guide_strand + '_' + presentallele + '_' + cas_obj.name + '_' + str(guide_len) + 'nt')
                    
            else:
                pass
    
    guides_df = pd.DataFrame({
        'chrom': chrom_list,
        'start': start_list,
        'end': end_list,
        'rsID': rsIDs,
        'SNP position': SNPcoordinates,
        'gRNA': guidelist,
        'PAM': PAMlist,
        'strand': strand,
        'variation': variation,
        'present allele': SNPallele,
        'guide ID': guidename,
        'SNP position in protospacer or PAM': pos_in_protospacer_or_pam
        
        })

    return guides_df


# COMBINE GUIDESDFs and ANNOTATE VARIANT INFO

"""
Input combined all_guides_df, and annotate in a new column whether each SNP is present in each vcf file input. Save alt AF frequency in a separate new column if SNP is in a population database vcf inputted.
"""
def all_guides_var_info(gensdict, guidesdf):
    
    guidesdf_copy = guidesdf.copy()
    new_guidesdf = guidesdf_copy.drop_duplicates()
    new_guidesdf = new_guidesdf.sort_values(by='start')

    for key in gensdict:
        in_vcf = new_guidesdf["SNP position"].isin(gensdict[key]["pos"])
        df_output = in_vcf.map({True: ("Y"), False: "N"})
        col_name = "In " + str(key)[8:]
        new_guidesdf[col_name] = list(df_output)

    alt_af_list = []

    for index, row in new_guidesdf.iterrows():
        snp_pos = row["SNP position"]
        af_values = []
        
        for vcf_name, vcf_df in gensdict.items():
            if "alt AF" in vcf_df.columns:
                matches = vcf_df[vcf_df["pos"] == snp_pos]
                if not matches.empty:
                    af_values.extend(matches["alt AF"].astype(str).tolist())
        
        alt_af_list.append(",".join(af_values) if af_values else "N/A")
    
    new_guidesdf["alt allele frequency"] = alt_af_list

    return new_guidesdf


# OFF-TARGET ANALYSIS FUNCTIONS

"""
3 functions to count exact matches in the genome for each guide in a guidesdf. Uses multi-processing for supposedly faster execution.
"""
genome = None
cas_obj = None

def init_worker(genome_path, cas_params):
    global genome, cas_obj
    genome = Fasta(genome_path)
    cas_obj = cas_params

def search_pattern(guide):
    total = 0
    total_rc = 0
    rc_guide = str(Seq(guide).reverse_complement())

    # Lookahead for overlapping matches, instead of overlapped=True. Probably can use overlapped=True too
    if cas_obj.is_three_prime:
        fwd_pattern = regex.compile(f"(?=({guide}{cas_obj.pam_three_prime}))", regex.IGNORECASE)
        rev_pattern = regex.compile(f"(?=({cas_obj.pam_five_prime}{rc_guide}))", regex.IGNORECASE)
    elif cas_obj.is_five_prime:
        fwd_pattern = regex.compile(f"(?=({cas_obj.pam_five_prime}{guide}))", regex.IGNORECASE)
        rev_pattern = regex.compile(f"(?=({rc_guide}{cas_obj.pam_three_prime}))", regex.IGNORECASE)
    else:
        raise ValueError("pam_orientation not defined")

    for chrom in genome.keys():
        seq = str(genome[chrom])
        total += len(fwd_pattern.findall(seq, overlapped=True))
        total_rc += len(rev_pattern.findall(seq, overlapped=True))

    return total + total_rc

def count_exact_matches(df, genome_fasta_path, cas_parameters, num_processes=None):
    guides = df['gRNA'].tolist()
    num_processes = num_processes or max(1, mp.cpu_count() - 1)

    # Start the multiprocessing pool
    with mp.Pool(
        num_processes,
        initializer=init_worker,
        initargs=(genome_fasta_path, cas_parameters)
    ) as pool:
        match_counts = pool.map(search_pattern, guides)

    df["exact matches in genome"] = match_counts
    return df


"""
Input: guidesdf and chromosome seq object. Counts 1 bp mismatches for each guide and stores in a new column in the dataframe.
"""
def one_mismatch(df, chseq):
    counts = []

    for idx in range(len(df)):
        guide = df['gRNA'].iloc[idx]
        guide_rc = str(Seq(guide).reverse_complement())

        chseq_str = str(chseq)

        #Generate mismatched guides. Only creates 1 bp mismatched guides, never searches for the exact guide. eg: [^T] means match any character except T 
        mismatch_guides = [guide[:j] + '[^' + guide[j] + ']' + guide[j+1:] for j in range(len(guide))] 
        mismatch_guides_rc = [guide_rc[:k] + '[^' + guide_rc[k] + ']' + guide_rc[k+1:] for k in range(len(guide_rc))] #create 1 bp mismatched guide set for reverse complemented guide

        count_per_guide = 0

        #Forward strand matches: guide + NGG
        for pat in mismatch_guides:
            regexp = regex.compile(pat + r'[ACGT]GG', regex.IGNORECASE)
            matches = regexp.findall(chseq_str, overlapped=True)
            count_per_guide += len(matches)

        #Reverse strand matches: CCN + rc_guide
        for pat_rc in mismatch_guides_rc:
            regexp_rc = regex.compile(r'CC[ACGT]' + pat_rc, regex.IGNORECASE)
            matches_rc = regexp_rc.findall(chseq_str, overlapped=True)
            count_per_guide += len(matches_rc)

        counts.append(count_per_guide)

    df['1 bp mismatches in ch1'] = counts
    return df



# SUMMARY STATS OF THE SINGLE-GUIDE LIBRARY

"""
Returns a list of variants that are targetable by your Cas enzyme (have PAM sites within 10 bp of the SNP), as well as the number of guides found for each SNP
"""
def targetable_vars(guidesdf):
    targetable_snps = guidesdf.loc[:, ['rsID', 'SNP position', 'alt allele frequency']]
    targetable_snps['no. of guides found (with ref or alt allele)'] = targetable_snps.groupby('SNP position')['SNP position'].transform('count')
    targetable_snps = targetable_snps.drop_duplicates(subset='SNP position', keep='first')
    targetable_snps = targetable_snps.reset_index() #makes new index column
    targetable_snps = targetable_snps.drop(['index'], axis=1)

    return targetable_snps



# PAIRING-RELATED FUNCTIONS

'''
Returns a list of two dataframes, split from the guides dataframe, by guides targeting  each allele, to enable pairing guides on the same allele.
'''
def split_phased(clean_guidesdf, phased_vcf, locus):
    phased_gensdf = create_gens(phased_vcf, locus, 'cell-line')
    
    #filtering out guides targeting SNPs that aren't heterozygous in the given vcf file
    filtered_guidesdf = clean_guidesdf[clean_guidesdf['SNP position'].isin(phased_gensdf['pos'])].copy()
    
    phasing_list = []

    for i in range(len(filtered_guidesdf)):

        # for the current row in guides df, getting the SNP allele present in the guide, and the position of the SNP
        allele = filtered_guidesdf['present allele'].iloc[i]
        pos = filtered_guidesdf['SNP position'].iloc[i]

        # extract row of VCF data of cell-line for that SNP
        match = phased_gensdf[phased_gensdf['pos'] == pos]

        #error handling if SNP is not found in VCF
        if match.empty:
            raise ValueError(
                f"Phased pairing requested but SNP at position {pos} in the guides dataframe is not present in the cell-line's VCF region."
            )

        # extracting the genotype from the extracted row, and splitting the genotype to create a list from the genotype data for easy parsing.
        gt = match['genotype'].iloc[0].split('|')

        # if the genotype wasn't able to be split by | (because the gt is unphased eg: 0/1), i.e., the pre-split and post-split value is the same, store its phasing information as "unphased". These guides will be included in both allele dataframes.
        if gt[0] == match['genotype'].iloc[0]:
            phasing_list.append('unphased')
        else:
            # storing the phasing of this version of the guide (ref or alt version of the SNP present) in this row of guidesdf
            if allele[3:6] == 'ref':
                if gt[0] == '0':
                    phasing_list.append('allele 1')
                elif gt[1] == '0':
                    phasing_list.append('allele 2')
                else:
                    phasing_list.append('error')
            else:
                if gt[0] == '1':
                    phasing_list.append('allele 1')
                elif gt[1] == '1':
                    phasing_list.append('allele 2')
                else:
                    phasing_list.append('error')

    filtered_guidesdf['phasing'] = phasing_list #adding phasing column to df

    #using this phasing, splitting the df into two guides dfs of phased guides.
    rows_allele1 = []
    rows_allele2 = []
    
    #for each row of the guidesdf
    for i in range(len(filtered_guidesdf)): 
        phasing = filtered_guidesdf['phasing'].iloc[i] #get phasing
        if phasing == 'allele 1':
            rows_allele1.append(filtered_guidesdf.iloc[i])
        elif phasing == 'allele 2':
            rows_allele2.append(filtered_guidesdf.iloc[i])
        else:
            rows_allele1.append(filtered_guidesdf.iloc[i])
            rows_allele2.append(filtered_guidesdf.iloc[i])


    #cleaning up each df
    filtered_guidesdf_allele1 = pd.DataFrame(rows_allele1)
    filtered_guidesdf_allele1 = filtered_guidesdf_allele1.reset_index()
    filtered_guidesdf_allele1 = filtered_guidesdf_allele1.drop(['index'], axis=1)
    
    filtered_guidesdf_allele2 = pd.DataFrame(rows_allele2)
    filtered_guidesdf_allele2 = filtered_guidesdf_allele2.reset_index()
    filtered_guidesdf_allele2 = filtered_guidesdf_allele2.drop(['index'], axis=1)
    
    return [filtered_guidesdf_allele1, filtered_guidesdf_allele2]


"""
Returns paired gRNA library with all possible guide pairs of a list of sgRNA
"""
def random_pair(guidesdf):
    data = []
    
    # iterate over unique combinations of gene sequences and IDs
    for i in range(len(guidesdf)):
        for j in range(i+1, len(guidesdf)):
            if guidesdf['SNP position'].iloc[i] != guidesdf['SNP position'].iloc[j]:
                data.append([guidesdf['SNP position'].iloc[i], guidesdf['alt allele frequency'].iloc[i], guidesdf['guide ID'].iloc[i], guidesdf['gRNA'].iloc[i], guidesdf['SNP position'].iloc[j], guidesdf['alt allele frequency'].iloc[j], guidesdf['guide ID'].iloc[j], guidesdf['gRNA'].iloc[j], (guidesdf['guide ID'].iloc[i] + '|' + guidesdf['guide ID'].iloc[j])])
            
    # convert to dataframe  
    pairings_df = pd.DataFrame(data, columns=['SNP position 1', 'alternate allele frequency 1', 'guide 1 ID', 'guide 1', 'SNP position 2', 'alternate allele frequency 2', 'guide 2 ID', 'guide 2', 'pair ID'])

    pairings_df = pairings_df.reset_index()
    pairings_df = pairings_df.drop(['index'], axis=1)

    #compute excision sizes for each guide pair and store in a new column
    excision_sizes = []
    
    for i in range(len(pairings_df)):
        guide_start1 = int((pairings_df['guide 1 ID'][i]).split("_")[1])
        guide_start2 = int((pairings_df['guide 2 ID'][i]).split("_")[1])
        excision_size = abs(guide_start1 - guide_start2)
        excision_sizes.append(excision_size)
        
    pairings_df['excision size'] = excision_sizes

    return pairings_df


"""
Returns paired gRNA library with guides paired to cause an excision about a fixed point. Takes as input a list of any number of such fixed points. Implemented to allow use cases such as excising certain exons or features of interest.
"""
def fixed_point_pair(guidesdf, points_list):
    data = []

    for point in points_list:
        point = int(point)
        for i in range(len(guidesdf)):
            if int(guidesdf['SNP position'].iloc[i]) <= point: #for snps below or equal to point
                for j in range(len(guidesdf)):
                    if int(guidesdf['SNP position'].iloc[j]) > point: #if next guide snp is above point
                        data.append([guidesdf['SNP position'].iloc[i], guidesdf['alt allele frequency'].iloc[i], guidesdf['guide ID'].iloc[i], guidesdf['gRNA'].iloc[i], guidesdf['SNP position'].iloc[j], guidesdf['alt allele frequency'].iloc[j], guidesdf['guide ID'].iloc[j], guidesdf['gRNA'].iloc[j], (guidesdf['guide ID'].iloc[i] + '|' + guidesdf['guide ID'].iloc[j])])
                    else:
                        'passing'
                        pass
    
            elif int(guidesdf['SNP position'].iloc[i]) >= point: #for snps above or equal to point
                for j in range(len(guidesdf)):
                    if int(guidesdf['SNP position'].iloc[j]) < point: #if next guide snp is below point
                        data.append([guidesdf['SNP position'].iloc[i], guidesdf['alt allele frequency'].iloc[i], guidesdf['guide ID'].iloc[i], guidesdf['gRNA'].iloc[i], guidesdf['SNP position'].iloc[j], guidesdf['alt allele frequency'].iloc[j], guidesdf['guide ID'].iloc[j], guidesdf['gRNA'].iloc[j], (guidesdf['guide ID'].iloc[i] + '|' + guidesdf['guide ID'].iloc[j])])
                    else:
                        'passing'
                        pass
                
        # convert to dataframe  
        fp_pairings_df = pd.DataFrame(data, columns=['SNP position 1', 'alternate allele frequency 1', 'guide 1 ID', 'guide 1', 'SNP position 2', 'alternate allele frequency 2', 'guide 2 ID', 'guide 2', 'pair ID'])

    fp_pairings_df = fp_pairings_df.drop_duplicates()
    fp_pairings_df = fp_pairings_df.reset_index()
    fp_pairings_df = fp_pairings_df.drop(['index'], axis=1)

    #compute excision sizes for each guide pair and store in a new column
    excision_sizes = []
    
    for i in range(len(fp_pairings_df)):
        guide_start1 = int((fp_pairings_df['guide 1 ID'][i]).split("_")[1])
        guide_start2 = int((fp_pairings_df['guide 2 ID'][i]).split("_")[1])
        excision_size = abs(guide_start1 - guide_start2)
        excision_sizes.append(excision_size)
    fp_pairings_df['excision size'] = excision_sizes

    return fp_pairings_df

def tiling_pair(guidesdf):
    guidesdf = guidesdf.sort_values(by=["SNP position","start"]).reset_index(drop=True)
    data = []

    for i in range(len(guidesdf)):
        snp_i = guidesdf['SNP position'].iloc[i]

        for j in range(i + 1, len(guidesdf)):
            snp_j = guidesdf['SNP position'].iloc[j]

            if snp_j == snp_i:
                continue  # same SNP, skip
            elif snp_j > snp_i:
                # First time SNP changes — pair i with all guides of the next SNP
                data.append([
                    snp_i,
                    guidesdf['alt allele frequency'].iloc[i],
                    guidesdf['guide ID'].iloc[i],
                    guidesdf['gRNA'].iloc[i],
                    snp_j,
                    guidesdf['alt allele frequency'].iloc[j],
                    guidesdf['guide ID'].iloc[j],
                    guidesdf['gRNA'].iloc[j],
                    guidesdf['guide ID'].iloc[i] + '|' + guidesdf['guide ID'].iloc[j]
                ])
            else:
                break  # if sorted, this won't happen — just for safety

            # now that we've reached the first different SNP, we continue checking if the next rows are also from that same SNP
            # and keep pairing i with them
            k = j + 1
            while k < len(guidesdf) and guidesdf['SNP position'].iloc[k] == snp_j:
                data.append([
                    snp_i,
                    guidesdf['alt allele frequency'].iloc[i],
                    guidesdf['guide ID'].iloc[i],
                    guidesdf['gRNA'].iloc[i],
                    guidesdf['SNP position'].iloc[k],
                    guidesdf['alt allele frequency'].iloc[k],
                    guidesdf['guide ID'].iloc[k],
                    guidesdf['gRNA'].iloc[k],
                    guidesdf['guide ID'].iloc[i] + '|' + guidesdf['guide ID'].iloc[k]
                ])
                k += 1
            break  # prevent pairing with SNPs beyond the next one

    # remove duplicates regardless of order (A|B == B|A)
    seen = set()
    unique_data = []
    for row in data:
        gid1, gid2 = row[2], row[6]
        pair_id = '|'.join(sorted([gid1, gid2]))
        if pair_id not in seen:
            seen.add(pair_id)
            row[8] = pair_id  # replace with ordered version
            unique_data.append(row)

    return pd.DataFrame(data, columns=[
        'SNP position 1', 'alternate allele frequency 1', 'guide 1 ID', 'guide 1',
        'SNP position 2', 'alternate allele frequency 2', 'guide 2 ID', 'guide 2',
        'pair ID'
    ])