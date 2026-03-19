# script to get the common variants in certain cpg sites
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from itertools import product

def main():

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--af_limit', type=float, required = True, help = "Allele frequency threshold.")
        parser.add_argument('--af_file_dir', type=str, required=True, help="Directory containing TGP_chr*_afs.txt files.")
        parser.add_argument('--exon_file', type=str, required = True, help = "Exon_information file.")
        #parser.add_argument('--acceptor_snp_region', type=str, required = True, help = "Window relative to exon start in which to look for snp. E.g., '3-21'.")
        parser.add_argument('--donor_snp_region', type=str, required = True, help="Window relative to exon start in which to look for snp. E.g., '4-21'.")
        parser.add_argument('--editing_window_size', type=str, required = True, help="Base editing window size. E.g., '4-8'.")
        parser.add_argument('--output_dir', type=str, required=True, help='Location for output files.')
        args = parser.parse_args()
        return args
    
    # parse args from input
    args=parse_args()
    af_limit=args.af_limit
    af_file_dir=args.af_file_dir
    exon_file=args.exon_file
    editing_window_size=args.editing_window_size
    donor_snp_region = args.donor_snp_region
    #acceptor_snp_region = args.acceptor_snp_region
    output_dir=args.output_dir

    # load exon file
    exon_df=pd.read_csv(exon_file, dtype={'chromsome_name':'str'}, index_col=0)
    # filter df to protein coding transcripts only
    protein_coding_biotypes=['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA']
    pc = exon_df[exon_df['transcript_biotype'].isin(protein_coding_biotypes)]

    # parse out regions that you'll search for snps in
    # acceptor_snp_region_lb = int(acceptor_snp_region.split('-')[0]) # lb for lower bound
    # acceptor_snp_region_ub = int(acceptor_snp_region.split('-')[1]) # ub for upper bound
    donor_snp_region_lb = int(donor_snp_region.split('-')[0])
    donor_snp_region_ub = int(donor_snp_region.split('-')[1])

    # find shared splice sites among each transcript in the df
    #acceptorSNP_ranges_byGene={}
    donorSNP_ranges_byGene={}
    for gene in pc.hgnc_symbol.unique():
        gene_df = pc[pc.hgnc_symbol==gene]
        # singleGene_transcript_acceptorSNP_ranges={}
        singleGene_transcript_donorSNP_ranges={}
        for transcript in gene_df.ensembl_transcript_id.unique():
            transcript_df=gene_df[gene_df.ensembl_transcript_id==transcript]
            cur_strand=transcript_df.strand.values[0]
            tdf_sorted = transcript_df.sort_values(by='rank', ascending=True)
            # get the splice acceptor and donor sites around each exon
            # cur_transcript_acceptor_ranges=[]
            cur_transcript_donor_ranges=[]
            for idx,row in tdf_sorted.iterrows():
                cur_exon=row.ensembl_exon_id
                if cur_strand==1:
                    exon_start = row.exon_chrom_start
                    exon_end = row.exon_chrom_end
                    # # snp range for splice acceptor alteration
                    # acceptor_range1_start = exon_start - acceptor_snp_region_ub
                    # acceptor_range1_end = exon_start - acceptor_snp_region_lb
                    # # second downstream range for acceptor (because of symmetry)
                    # acceptor_range2_start = exon_start + acceptor_snp_region_lb
                    # acceptor_range2_end = exon_start + acceptor_snp_region_ub
                    # snp range for splice donor alteration
                    donor_range_start = exon_end - donor_snp_region_ub
                    donor_range_end = exon_end - donor_snp_region_lb
                else: # if current_strand = -1
                    exon_start = row.exon_chrom_end
                    exon_end = row.exon_chrom_start
                    # # snp range for splice acceptor alteration
                    # acceptor_range1_start = exon_start - acceptor_snp_region_ub
                    # acceptor_range1_end = exon_start - acceptor_snp_region_lb
                    # # second downstream range for acceptor
                    # acceptor_range2_start = exon_start + acceptor_snp_region_lb
                    # acceptor_range2_end = exon_start + acceptor_snp_region_ub
                    # snp range for splice donor alteration
                    donor_range_start = exon_end + donor_snp_region_lb
                    donor_range_end = exon_end + donor_snp_region_ub

                # great, now store these ranges
                # cur_transcript_acceptor_ranges.append([cur_exon, acceptor_range1_start, acceptor_range1_end])
                # cur_transcript_acceptor_ranges.append([cur_exon, acceptor_range2_start, acceptor_range2_end])
                cur_transcript_donor_ranges.append([cur_exon, donor_range_start, donor_range_end])

            # store the current acceptor/donor ranges for the 
            # singleGene_transcript_acceptorSNP_ranges[transcript] = pd.DataFrame(cur_transcript_acceptor_ranges, columns=['ensembl_exon_id','acceptor_range_start', 'acceptor_range_end'])
            singleGene_transcript_donorSNP_ranges[transcript] = pd.DataFrame(cur_transcript_donor_ranges, columns=['ensembl_exon_id','donor_range_start', 'donor_range_end'])

        # finally, store these for the gene overall
        # acceptorSNP_ranges_byGene[gene] = singleGene_transcript_acceptorSNP_ranges
        donorSNP_ranges_byGene[gene] = singleGene_transcript_donorSNP_ranges


    # Now, get the regions that are common among these
    # for donor ranges
    universal_donor_snp_regions_byGene={}
    for gene, transcript_dict in donorSNP_ranges_byGene.items():
        # get number of transcripts for the current gene
        transcript_list = list(transcript_dict.keys())
        num_transcripts = len(transcript_list)
        
        # initialize the first range to compare to as the first transcript's donor regions
        if len(transcript_list)>0:
            transcript1 = transcript_list[0]
            t1_df = transcript_dict[transcript_list[0]]
            t1_df['start_end'] = t1_df[['donor_range_start','donor_range_end']].apply(tuple,axis=1)
            t1_donor_loc_list = list(t1_df['start_end'])

            # have an if statement here checking if there is only one transcript for the gene (otherwise don't need to run all of the stuff below)
            for transcript2, t2_df in transcript_dict.items():
                if len(t1_donor_loc_list)!=0:
                    if transcript1!=transcript2:
                        # if the transcripts we're comparing aren't the same, compare their donor ranges, one exon of t1 at a time
                        t2_copy = t2_df.copy()
                        t2_num_ranges = len(t2_copy.ensembl_exon_id)
                        
                        # turn start and end into tuples
                        t2_copy['start_end'] = t2_copy[['donor_range_start','donor_range_end']].apply(tuple,axis=1)
                        t2_donor_loc_list = list(t2_copy['start_end'])
                        
                        # create a product to get every combination of exon ranges that can be compared
                        range_combos = list(product(t2_donor_loc_list, t1_donor_loc_list))
                        combos_df = pd.DataFrame(range_combos,columns=['t2_vals','t1_vals'])
                        combos_df[['t2a', 't2b']] = combos_df.t2_vals.tolist()
                        combos_df[['t1a','t1b']] = combos_df.t1_vals.tolist()
                        combos_df = combos_df[['t2a','t2b','t1a','t1b']]
                        
                        # calculate overlap between each combo
                        combos_df['overlap_start']=combos_df[['t1a','t2a']].max(axis=1)
                        combos_df['overlap_end']=combos_df[['t1b','t2b']].min(axis=1)
                        combos_df['overlap'] = combos_df.overlap_end - combos_df.overlap_start
                        combos_df['overlap_bool'] = combos_df.overlap>0 # we'll consider 0 overlapping because this means a single bp is shared
                        
                        # now let's filter to just the regions that overlap
                        o_regions = combos_df[combos_df.overlap_bool==True]
                        # convert overlapping regions back into tuples
                        o_regions = o_regions.assign(overlapping_regions = o_regions[['overlap_start','overlap_end']].apply(tuple,axis=1))
        
                        # now reset t1_exon_loc_list to compare to the rest of the transcripts
                        t1_donor_loc_list = list(o_regions['overlapping_regions'])
                else:
                    break
        else:
            t1_donor_loc_list=[]

        # save the 'universal donor snp regions' that we found
        universal_donor_snp_regions_byGene[gene] = t1_donor_loc_list
        with open(output_dir + '/ubiq_regions/ubiq_donorRegions_ALL_chroms.pkl','wb') as file:
            pickle.dump(universal_donor_snp_regions_byGene, file)


    ## -----------------------------------------------------------------
    # Now, find the common variants in these regions for each gene
    
    # load the vcf files for biallelic snps
    vcf_dict={}
    chroms = pc.chromosome_name.unique()
    for chrom in chroms:
        af_filename = os.path.join(af_file_dir, 'TGP_chr' + str(chrom) + '_afs.txt')
        cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
        cur_chrom_TGP_afs = cur_chrom_TGP_afs[(cur_chrom_TGP_afs.af>=af_limit) & (cur_chrom_TGP_afs.af<=1-af_limit)] # filter to af threshold
        vcf_dict[chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt', 'af']]

    # loop over the shared regions to check which ones have a snp in them at the specified allele frequency
    num_common_vars_in_donor_regions=[]
    genes=[]
    chroms=[]
    common_var_donor_info = {}
    for gene, region_list in universal_donor_snp_regions_byGene.items():
        genes.append(gene)
        cur_chrom = pc[pc.hgnc_symbol==gene]['chromosome_name'].values[0]
        chroms.append(cur_chrom)
        cur_vcf = vcf_dict[cur_chrom]
        cur_gene_common_var_donor_info=[]
        if len(region_list)>0:
            # turn list into a data frame
            for region in region_list:
                region_start=region[0]
                region_end=region[1]
                in_donor_regions=[x for x in list(cur_vcf.pos) if (x>=region_start and x<=region_end)]
                # get the allele frequency information for these positions too
                pos_afs = cur_vcf[cur_vcf['pos'].isin(in_donor_regions)][['pos','af']]
                for idx,row in pos_afs.iterrows():
                    cur_gene_common_var_donor_info.append([row.pos,row.af])
            num_common_vars_in_donor_regions.append(len(cur_gene_common_var_donor_info))
        else:
            num_common_vars_in_donor_regions.append(None)
            cur_gene_common_var_donor_info = ['No shared donor regions']
        common_var_donor_info[gene] = cur_gene_common_var_donor_info
    # combine summary information into data frame
    donor_cv_df = pd.DataFrame({
        'gene':genes,
        'num_common_vars_in_donor_regions':num_common_vars_in_donor_regions,
        'chrom':chroms
    })
    # save the information here
    donor_cv_df.to_csv(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_summary.txt" ,sep='\t')
    donor_cv_df.to_csv(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" ,sep='\t', index=False, header=False)
    with open(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl",'wb') as file:
        pickle.dump(common_var_donor_info,file)




# --------------------------------
if __name__ == '__main__':
    main()
