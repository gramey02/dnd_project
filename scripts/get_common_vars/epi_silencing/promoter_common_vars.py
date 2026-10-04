import argparse
import os
import sys
import math
import numpy as np
import pandas as pd
import pickle
import pysam

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def atomic_pickle_dump(obj, path):
    # write-then-rename so concurrent array tasks never see a partially-written file under its final name
    tmp_path = path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump(obj, f)
    os.rename(tmp_path, path)

def atomic_to_csv(df, path, **kwargs):
    tmp_path = path + '.tmp'
    df.to_csv(tmp_path, **kwargs)
    os.rename(tmp_path, path)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def calc_gc_content(chrom = None,
                    start_region=None,
                    end_region=None,
                    fasta_open_obj=None):
    """
    This function calculates the gc content for a particular sequence of DNA.
    Parameters:
        chrom: str
            Chromosome that the sequence falls in
        start_region: int
            start of region (e.g., start of a promoter)
        end_region: int
            end of region (e.g., end of promoter)
        fasta_open_obj: pysam object
            object created with the command "fasta_open = pysam.Fastafile(fasta_filepath)", where fasta_filepath is the path to the human reference genome
    Returns:
        gc content (percentage of total region)
    """
    # get seq
    seq = fasta_open_obj.fetch(start = start_region, end = end_region, region = str(chrom))
    # find out how many Cs and Gs are in the seq string
    letter_to_find="c"
    c_count = seq.lower().count(letter_to_find.lower())
    letter_to_find="g"
    g_count = seq.lower().count(letter_to_find.lower())
    total_gc = c_count + g_count
    gc_content_percent = (total_gc/(end_region-start_region))*100
    return gc_content_percent

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_file', type = str, required = True, help = 'CpG locations filepath.')
    parser.add_argument('--gt_file', type = str, required = True, help = 'Filepath that contains gene, transcript, and exon information for each gene.')
    parser.add_argument('--intersection', type = int, required = True, help = 'Amount of base pairs to consider two regions overlapping.')
    parser.add_argument('--promoter_ud', type = int, required = True, help = 'Distance upstream of transcription start site that promoter begins at.')
    parser.add_argument('--promoter_dd', type = int, required = True, help = 'Distance downstream of transcription start site that promoter begins at.')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory to save information to.')
    parser.add_argument('--af_limit', type=float, required=True, help="Allele frequency limit.")
    parser.add_argument('--af_file_dir', type=str, required=True, help="Directory containing TGP_chr*_afs.txt files.")
    parser.add_argument('--gc_threshold', type=float, required=False, help='GC content threshold to filter promoters by.')
    parser.add_argument('--use_islands', type=str,required=True, help='Boolean that tells the script whether to filter promoter regions by overlap with a CpG island or by average GC threshold of the promoter region')
    parser.add_argument('--ref_genome_fasta', type=str, required=True, help='Reference genome fasta file path.')
    parser.add_argument('--chrom', type=str, required=True, help='Chromosome to process.')
    parser.add_argument('--total_num_chroms', type=int, required=True, help='Total number of unique chromosomes for gene set of interest.')
    args = parser.parse_args()
    return args

def main():
    # get the passed args
    args = parse_args()
    cpg_file = args.cpg_file
    gt_file = args.gt_file
    acceptable_intersection = args.intersection
    promoter_ud = args.promoter_ud
    promoter_dd = args.promoter_dd
    output_dir=args.output_dir
    af_limit=args.af_limit
    af_file_dir=args.af_file_dir
    gc_thresh=args.gc_threshold
    use_islands=args.use_islands
    ref_genome_fasta=args.ref_genome_fasta
    chrom=args.chrom
    total_num_chroms=args.total_num_chroms

    # load cpgs
    cpg_df = pd.read_table(cpg_file)
    cpg_df.columns=['chrom', 'chromStart', 'chromEnd', 'cpg_name', 'length', 'cpgNum', 'gcNum', 'perCpg', 'perGc', 'obsExp']
    chroms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    cpg_df=cpg_df[cpg_df['chrom'].isin(chroms)]
    # only need this chromosome's cpg islands
    cpg_df=cpg_df[cpg_df['chrom']=='chr'+str(chrom)]

    # load gene & transcript information
    dtype_dict_small={'chromosome_name':'str'} # set type of one of the columns, since it's ambiguous atm
    gt_df = pd.read_csv(gt_file,index_col=0,dtype=dtype_dict_small)
    gt_df=(gt_df.drop(labels=["gene_biotype","ensembl_peptide_id","exon_chrom_start","exon_chrom_end","ensembl_exon_id","is_constitutive","rank", "cds_start","cds_end", "genomic_coding_start", "genomic_coding_end"],axis=1)).drop_duplicates()
    # the input file to this script may already represent transcripts filtered by expression (those with more rare expression should've been filtered out)
    # let's also filter to protein coding transcripts only
    protein_coding_biotypes=['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA']
    gt_df=gt_df[gt_df.transcript_biotype.isin(protein_coding_biotypes)].reset_index(drop=True) # ge_cds = gene_exons, coding sequences
    # only process this chromosome's genes (this script is run as one SGE array task per chromosome)
    gt_df=gt_df[gt_df.chromosome_name==chrom].reset_index(drop=True)

    # get range of promoter for each transcript
    promoter_starts=[]
    promoter_ends=[]
    for idx,row in gt_df.iterrows():
        cur_transcript=row.ensembl_transcript_id
        # get promoter bounds for the transcript
        if row.strand==1:
            promoter_start = row.transcription_start_site - promoter_ud
            promoter_end = row.transcription_start_site + promoter_dd
        else:
            promoter_start = row.transcription_start_site - promoter_dd
            promoter_end = row.transcription_start_site + promoter_ud

        promoter_starts.append(promoter_start)
        promoter_ends.append(promoter_end)
    gt_df['promoter_start'] = promoter_starts
    gt_df['promoter_end'] = promoter_ends

    # ---------------------------------------------
    # calculate the CpG content of the promoter region of a gene
    fasta_filepath = ref_genome_fasta # "/wynton/group/capra/data/hg38_fasta/2022-03-14/hg38.fa.gz"
    fasta_open = pysam.Fastafile(fasta_filepath) # open the fasta file for an individual sample
    gc_content=[]
    for idx,row in gt_df.iterrows():
        gc_content.append(calc_gc_content(chrom=row.chromosome_name,start_region=row.promoter_start, end_region=row.promoter_end, fasta_open_obj=fasta_open))
    gt_df['promoter_GCcontent_percent'] = gc_content

    # average gc content across transcript promoters per gene
    grouped = pd.DataFrame(gt_df.groupby('ensembl_gene_id')['promoter_GCcontent_percent'].mean())
    grouped.rename(columns={'promoter_GCcontent_percent':'avg_promoter_GC_percent'},inplace=True)
    gt_df = gt_df.merge(grouped,on='ensembl_gene_id',how='left')
    # -----------------------------------------------
    # now get the shared promoter regions for each gene
    gene_list = []
    overlap_len_list = []
    overlap_bool = []
    overlap_start_list = []
    overlap_end_list = []
    overlap_end = []
    for gene in gt_df.ensembl_gene_id.unique():
        gene_list.append(gene)
        gene_df = gt_df[gt_df.ensembl_gene_id==gene]
        overlap_start = max(gene_df.promoter_start)
        overlap_end = min(gene_df.promoter_end)
        overlap_len = overlap_end-overlap_start
        if overlap_len>=0:
            overlap_bool.append(1)
            overlap_start_list.append(overlap_start)
            overlap_end_list.append(overlap_end)
        else:
            overlap_bool.append(0)
            overlap_start_list.append(None)
            overlap_end_list.append(None)
        overlap_len_list.append(overlap_len)
    overlap_data = pd.DataFrame({
        'ensembl_gene_id':gene_list,
        'transcripts_have_overlapping_promoters':overlap_bool,
        'overlap_length':overlap_len_list,
        'overlap_start':overlap_start_list,
        'overlap_end':overlap_end_list
    })
    gt_df = gt_df.merge(overlap_data,on='ensembl_gene_id',how='left')
    #---------------------------------------------------------
    # now let's check, for each promoter if it overlaps with a cpg site
    gene_cpg_dict={}
    overlaps_cpg_island=[]
    overlap_amount=[]
    gene_list=[]
    # filter to just the genes that have shared promoter regions
    promoter_sharing_genes = gt_df[gt_df['transcripts_have_overlapping_promoters']==1]
    promoter_sharing_genes = (promoter_sharing_genes.drop(labels=['ensembl_transcript_id', 'hgnc_symbol','start_position', 'end_position','transcript_start', 'transcript_end', 'transcription_start_site',
        'transcript_length', 'transcript_is_canonical', 'transcript_tsl',
        'transcript_biotype', 'promoter_start', 'promoter_end',
        'promoter_GCcontent_percent', 'avg_promoter_GC_percent',
        'transcripts_have_overlapping_promoters', 'overlap_length'],axis=1)).drop_duplicates()
    promoter_sharing_genes.rename(columns={'overlap_start':'promoter_start','overlap_end':'promoter_end'},inplace=True)
    #counter=0
    for idx,row in promoter_sharing_genes.iterrows():
        gene_list.append(row.ensembl_gene_id)
        # filter to chromosome of the current transcript
        chrom_cpgs=cpg_df[cpg_df.chrom=='chr'+str(row.chromosome_name)]
        # add promoter information to the cpg df so you can easily calculate values for each cpg
        chrom_cpgs = chrom_cpgs.assign(promoter_start=([row.promoter_start]*len(chrom_cpgs.cpg_name)))
        chrom_cpgs = chrom_cpgs.assign(promoter_end=([row.promoter_end]*len(chrom_cpgs.cpg_name)))
        # calculate the overlap start and end
        chrom_cpgs = chrom_cpgs.assign(overlap_start=(chrom_cpgs[['promoter_start','chromStart']].max(axis=1)))
        chrom_cpgs = chrom_cpgs.assign(overlap_end=(chrom_cpgs[['promoter_end','chromEnd']].min(axis=1)))

        # calculate overlap
        chrom_cpgs = chrom_cpgs.assign(overlap=(chrom_cpgs.overlap_end-chrom_cpgs.overlap_start))
        # check where the overlap is less than the acceptable amount
        chrom_cpgs = chrom_cpgs.assign(overlap_bool=(chrom_cpgs.overlap>=acceptable_intersection))
        # filter to get only the CpGs that overlap
        overlapping_cpgs = chrom_cpgs[chrom_cpgs.overlap_bool==True]
        # save these for the gene

        if len(overlapping_cpgs.cpg_name)==0:
            gene_cpg_dict[row.ensembl_gene_id]=None
            overlaps_cpg_island.append(0)
            overlap_amount.append(None)
        else:
            if len(overlapping_cpgs.cpg_name)>1:
                gene_cpg_dict[row.ensembl_gene_id]=[list(overlapping_cpgs.cpg_name),list(overlapping_cpgs.chromStart),list(overlapping_cpgs.chromEnd)]
                overlaps_cpg_island.append(1)
                overlap_amount.append(np.mean(overlapping_cpgs.overlap))
            else:
                gene_cpg_dict[row.ensembl_gene_id]=[overlapping_cpgs.cpg_name.values[0], overlapping_cpgs.chromStart.values[0], overlapping_cpgs.chromEnd.values[0]]
                overlaps_cpg_island.append(1)
                overlap_amount.append(overlapping_cpgs.overlap.values[0])
        #counter+=1
    ubiq_dir = output_dir + '/ubiq_regions/'
    atomic_pickle_dump(gene_cpg_dict, ubiq_dir + 'dhs_CpGIsland_Overlap_chr' + chrom + '_dict.pkl')
    # merging across chromosomes now happens in a separate post-array step (merge_promoter_common_vars.py)
    # merge with larger df
    cpg_data = pd.DataFrame({
        'ensembl_gene_id':gene_list,
        'overlaps_CpGIsland':overlaps_cpg_island,
        'overlap_amount':overlap_amount
    })
    gt_df = gt_df.merge(cpg_data,on='ensembl_gene_id',how='left')
    #--------------------------------------------------------------------
    # now get things on a gene level by dropping columns and removing duplicates
    gt_df_filt = gt_df.drop(labels=['ensembl_transcript_id','promoter_start','promoter_end','start_position', 'end_position', 'strand',
        'transcript_start', 'transcript_end', 'transcription_start_site',
        'transcript_length', 'transcript_is_canonical', 'transcript_tsl',
        'transcript_biotype', 'promoter_GCcontent_percent'],axis=1).drop_duplicates()
    # rename columns for clarity
    gt_df_filt.rename(columns={'overlap_length':'shared_promoter_overlap_length',
                    'overlap_start':'shared_promoter_start',
                    'overlap_end':'shared_promoter_end',
                    'overlap_amount':'CpGIsland_overlap_amount'},inplace=True)
    #----------------------------------------------------------------------
    # convert the form of the promoter regions
    shared_promoter_regions={}
    gene_chrom_map={}
    #genes_w_shared_promoter_regions=[]
    for idx,row in gt_df_filt.iterrows():
        gene_chrom_map[row.hgnc_symbol]=row.chromosome_name
        if row.transcripts_have_overlapping_promoters==1:
            shared_promoter_regions[row.hgnc_symbol]=[tuple((row.shared_promoter_start, row.shared_promoter_end))]
            #genes_w_shared_promoter_regions.append(row.gene)
        else:
            shared_promoter_regions[row.hgnc_symbol]=['No shared promoter regions']
    # save the shared promoter region information
    atomic_pickle_dump(shared_promoter_regions, ubiq_dir + "ubiq_promoters_chr" + chrom + ".pkl")
    # merging across chromosomes now happens in a separate post-array step (merge_promoter_common_vars.py)
    #------------------------------------------------------------------------
    # load allele frequency files for common variant localization
    vcf_dict={}
    genes_chroms = list(set(gene_chrom_map.values()))
    for gene_chrom in genes_chroms:
        af_filename = os.path.join(af_file_dir, 'TGP_chr' + gene_chrom + '_afs.txt')
        cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
        cur_chrom_TGP_afs=cur_chrom_TGP_afs[(cur_chrom_TGP_afs.af>=af_limit) & (cur_chrom_TGP_afs.af<=1-af_limit)]
        cur_chrom_TGP_afs.drop_duplicates(subset='pos',keep=False,inplace=True)
        # get a list of common variant positions in and around the gene
        vcf_dict[gene_chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt', 'af']]
    #------------------------------------------------------------------------
    # finally, locate the common vars in the shared promoter regions
    common_var_info={}
    num_common_vars_in_promoters=[]
    gene_list=[]
    for gene,promoter_regions in shared_promoter_regions.items():
        gene_list.append(gene)
        cur_chrom=gene_chrom_map[gene]
        af_pass=vcf_dict[cur_chrom]
        pos=list(af_pass.pos)
        cur_gene_common_var_info=[]
        if len(promoter_regions)!=0:
            if type(promoter_regions[0])!=type('example_string'):
                # turn list into a data frame
                for region in promoter_regions:
                    region_start=region[0]
                    region_end=region[1]
                    in_promoter_regions=[x for x in pos if (x>=region_start and x<=region_end)]
                    # get the allele frequency information for these positions too
                    pos_afs = af_pass[af_pass['pos'].isin(in_promoter_regions)][['pos','af']]
                    for idx,row in pos_afs.iterrows():
                        cur_gene_common_var_info.append([row.pos,row.af])
                num_common_vars_in_promoters.append(len(cur_gene_common_var_info))
            else:
                num_common_vars_in_promoters.append(0)
                cur_gene_common_var_info = ['No shared promoter regions']
        else:
            num_common_vars_in_promoters.append(0)
            cur_gene_common_var_info = ['Promoter length = 0']
        if len(cur_gene_common_var_info)==0:
            cur_gene_common_var_info=['No variants pass MAF threshold']
        common_var_info[gene] = cur_gene_common_var_info
    #----------------------------------------------------------------------
    # filter either by CpG Island overlap or average promoter GC content
    final_var_info={}
    gene_list=[]
    num_common_vars_in_promoters=[]
    if (use_islands=='False') or (use_islands==False) or (use_islands==0):
        for gene,var_list in common_var_info.items():
            gene_list.append(gene)
            cur_gene_gc_info=gt_df_filt[gt_df_filt.hgnc_symbol==gene]
            if (var_list==['No shared promoter regions'] or var_list==['No variants pass MAF threshold'] or var_list==['Promoter length = 0']):
                final_var_info[gene]=var_list
                num_common_vars_in_promoters.append(0)
            else:
                if (cur_gene_gc_info.avg_promoter_GC_percent.values[0])>=gc_thresh:
                    final_var_info[gene]=var_list
                    num_common_vars_in_promoters.append(len(var_list))
                else:
                    final_var_info[gene]=['Does not pass GC content filter']
                    num_common_vars_in_promoters.append(0)
    else:
        for gene,var_list in common_var_info.items():
            gene_list.append(gene)
            cur_gene_gc_info=gt_df_filt[gt_df_filt.hgnc_symbol==gene]
            if (var_list==['No shared promoter regions'] or var_list==['No variants pass MAF threshold'] or var_list==['Promoter length = 0']):
                final_var_info[gene]=var_list
                num_common_vars_in_promoters.append(0)
            else:
                if math.isnan(cur_gene_gc_info.overlaps_CpGIsland.values[0])==False:
                    if (cur_gene_gc_info.overlaps_CpGIsland.values[0])==1:
                        final_var_info[gene]=var_list
                        num_common_vars_in_promoters.append(len(var_list))
                    else:
                        final_var_info[gene]=['Does not overlap CpG island']
                        num_common_vars_in_promoters.append(0)
                else:
                    final_var_info[gene]=['Does not overlap CpG island']
                    num_common_vars_in_promoters.append(0)

    #----------------------------------------------------------------------
    # save the information
    # combine summary information into data frame
    promoter_cv_df = pd.DataFrame({
        'gene':gene_list,
        'num_common_vars_in_shared_promoters':num_common_vars_in_promoters
    })
    # reformat
    promoter_cv_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
    promoter_df_merged = promoter_cv_df.merge(gt_df_filt[['hgnc_symbol','ensembl_gene_id','chromosome_name','transcripts_have_overlapping_promoters','avg_promoter_GC_percent','shared_promoter_overlap_length','overlaps_CpGIsland']], on='hgnc_symbol',how='left')
    # save this chromosome's results
    common_vars_savepath = output_dir + '/ubiq_region_CommonVars/'
    atomic_to_csv(promoter_df_merged, common_vars_savepath + 'CommonVars_chr' + chrom + '_summary.txt', sep='\t')
    atomic_pickle_dump(final_var_info, common_vars_savepath + 'CommonVars_chr' + chrom + '_dict.pkl')
    # merging across chromosomes now happens in a separate post-array step (merge_promoter_common_vars.py)


# -----------------------------------------------------------------------------------------------
if __name__=='__main__':
    main()
