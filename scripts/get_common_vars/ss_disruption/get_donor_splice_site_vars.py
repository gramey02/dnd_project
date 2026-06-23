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
        parser.add_argument('--donor_snp_region', type=str, required = True, help="Window relative to exon start in which to look for snp. E.g., '4-21'.")
        parser.add_argument('--editing_window_size', type=str, required = True, help="Base editing window size. E.g., '4-8'.")
        parser.add_argument('--output_dir', type=str, required=True, help='Location for output files.')
        parser.add_argument('--base_editor', type=str, required=False, help='Base editor to filter variants by.')
        args = parser.parse_args()
        return args
    
    # -----------------------------
    # helpers
    # -----------------------------
    def get_snps_in_range(vcf_df, start, end):
        return vcf_df[(vcf_df.pos >= start) & (vcf_df.pos <= end)].copy()


    def compute_distance(pos, boundary, strand):
        # signed distance relative to exon boundary
        if strand == 1:
            return pos - boundary
        else:
            return boundary - pos
        
    def snp_in_regions(snp_pos, regions):
        return any(start <= snp_pos <= end for start, end in regions)
    
    def explode_parallel_lists(df):
        out_rows = []

        for _, row in df.iterrows():

            snps = row["snps"]
            afs = row["snps_af"]
            dists = row["snp_distances"]

            # skip empty entries safely
            if len(snps) == 0:
                continue

            for snp, af, dist in zip(snps, afs, dists):
                new_row = row.to_dict()

                new_row["snp"] = snp
                new_row["snp_af"] = af
                new_row["snp_distance"] = dist

                # optionally drop original list columns
                del new_row["snps"]
                del new_row["snps_af"]
                del new_row["snp_distances"]

                out_rows.append(new_row)

        return pd.DataFrame(out_rows)
    
    # ------------------------
    # main start
    # ------------------------
    
    # parse args from input
    args=parse_args()
    af_limit=args.af_limit
    af_file_dir=args.af_file_dir
    exon_file=args.exon_file
    donor_snp_region = args.donor_snp_region
    output_dir=args.output_dir
    base_editor=args.base_editor

    # -----------------------------
    # load exons
    # -----------------------------
    exon_df = pd.read_csv(exon_file, index_col=0)
    pc = exon_df[exon_df['transcript_biotype'].isin(
        ['protein_coding','nonsense_mediated_decay','non_stop_decay','lncRNA','miRNA']
    )]
    # separate out lower and upper bound coordinates to look for snps in
    lb, ub = map(int, donor_snp_region.split('-'))
    # specific base editor range assignment
    ABE_LB, ABE_UB = lb, ub-1
    CBE_LB, CBE_UB = lb+1, ub

    # -----------------------------
    # load SNPs per chromosome
    # -----------------------------
    vcf_dict = {}
    for chrom in pc.chromosome_name.unique():
        f = os.path.join(af_file_dir, f"TGP_chr{chrom}_afs.txt")
        df = pd.read_csv(
            f, sep=' ',
            names=['chrom','pos','ref','alt','ac','an','af','afr_af','amr_af','eas_af','eur_af','sas_af']
        )

        df = df[(df.af >= af_limit) & (df.af <= 1 - af_limit)]
        vcf_dict[chrom] = df[['chrom','pos','af']]

    # =========================================================
    # STEP 1: build exon-level annotated objects (DONOR VERSION)
    # =========================================================

    gene_objects = {}

    for gene in pc.hgnc_symbol.unique():

        gene_df = pc[pc.hgnc_symbol == gene]
        gene_records = []

        for transcript in gene_df.ensembl_transcript_id.unique():

            tdf = gene_df[gene_df.ensembl_transcript_id == transcript]
            tdf = tdf.sort_values("rank")
            strand=tdf.strand.values[0]

            for _, row in tdf.iterrows():

                exon_id = row.ensembl_exon_id

                # donor window (single direction logic; strand irrelevant here)
                if strand==1:
                    exon_start = row.exon_chrom_start
                    exon_end = row.exon_chrom_end
                    # snp range for splice donor alteration
                    donor_start = exon_end - ub
                    donor_end = exon_end - lb

                    # ABE window
                    gene_records.append({
                        "gene": gene,
                        "transcript": transcript,
                        "exon_id": exon_id,
                        "strand": strand,
                        "range_start": exon_end - ABE_UB,
                        "range_end": exon_end - ABE_LB,
                        "editor": "ABE",
                        "exon_boundary": exon_end,
                        "boundary_type": "donor"
                    })

                    # CBE window
                    gene_records.append({
                        "gene": gene,
                        "transcript": transcript,
                        "exon_id": exon_id,
                        "strand": strand,
                        "range_start": exon_end - CBE_UB,
                        "range_end": exon_end - CBE_LB,
                        "editor": "CBE",
                        "exon_boundary": exon_end,
                        "boundary_type": "donor"
                    })


                else: # if current_strand = -1
                    exon_start = row.exon_chrom_end
                    exon_end = row.exon_chrom_start
                    # snp range for splice donor alteration
                    donor_start = exon_end + lb
                    donor_end = exon_end + ub

                    # ABE window
                    gene_records.append({
                        "gene": gene,
                        "transcript": transcript,
                        "exon_id": exon_id,
                        "strand": strand,
                        "range_start": exon_end + ABE_LB,
                        "range_end": exon_end + ABE_UB,
                        "editor": "ABE",
                        "exon_boundary": exon_end,
                        "boundary_type": "donor"
                    })

                    # CBE window
                    gene_records.append({
                        "gene": gene,
                        "transcript": transcript,
                        "exon_id": exon_id,
                        "strand": strand,
                        "range_start": exon_end + CBE_LB,
                        "range_end": exon_end + CBE_UB,
                        "editor": "CBE",
                        "exon_boundary": exon_end,
                        "boundary_type": "donor"
                    })

        gene_objects[gene] = pd.DataFrame(gene_records)

    # =========================================================
    # STEP 2: attach SNPs + distances
    # =========================================================
    enriched_gene_objects = {}

    for gene, df in gene_objects.items():

        chrom = pc[pc.hgnc_symbol == gene]['chromosome_name'].values[0]
        vcf = vcf_dict[chrom]

        enriched_rows = []

        for _, r in df.iterrows():

            snps = get_snps_in_range(vcf, r.range_start, r.range_end)

            snp_list = snps.pos.tolist()
            af_list = snps.af.tolist()

            distances = [
                compute_distance(p, r.exon_boundary, r.strand)
                for p in snp_list
            ]

            enriched_rows.append({
                **r.to_dict(),
                "snps": snp_list,
                "snps_af": af_list,
                "snp_distances": distances
            })

        enriched_gene_objects[gene] = pd.DataFrame(enriched_rows)
    
    # =========================================================
    # STEP 3: optional base editor filtering
    # =========================================================
    final_gene_objects = {}
    for gene, df in enriched_gene_objects.items():

        if base_editor in ["ABE", "CBE"]:
            df = df[df.editor == base_editor]

        final_gene_objects[gene] = df

    
    # =========================================================
    # STEP 4: RUN ORIGINAL CODE to get shared donor snp ranges
    # =========================================================
    # find shared splice sites among each transcript in the df
    donor_snp_region_lb = int(donor_snp_region.split('-')[0])
    donor_snp_region_ub = int(donor_snp_region.split('-')[1])
    donorSNP_ranges_byGene={}
    for gene in pc.hgnc_symbol.unique():
        gene_df = pc[pc.hgnc_symbol==gene]
        singleGene_transcript_donorSNP_ranges={}
        for transcript in gene_df.ensembl_transcript_id.unique():
            transcript_df=gene_df[gene_df.ensembl_transcript_id==transcript]
            cur_strand=transcript_df.strand.values[0]
            tdf_sorted = transcript_df.sort_values(by='rank', ascending=True)
            # get the splice acceptor and donor sites around each exon
            cur_transcript_donor_ranges=[]
            for idx,row in tdf_sorted.iterrows():
                cur_exon=row.ensembl_exon_id
                if cur_strand==1:
                    exon_start = row.exon_chrom_start
                    exon_end = row.exon_chrom_end
                    # snp range for splice donor alteration
                    donor_range_start = exon_end - donor_snp_region_ub
                    donor_range_end = exon_end - donor_snp_region_lb
                else: # if current_strand = -1
                    exon_start = row.exon_chrom_end
                    exon_end = row.exon_chrom_start
                    # snp range for splice donor alteration
                    donor_range_start = exon_end + donor_snp_region_lb
                    donor_range_end = exon_end + donor_snp_region_ub

                # great, now store these ranges
                cur_transcript_donor_ranges.append([cur_exon, donor_range_start, donor_range_end])

            # store the current acceptor/donor ranges for the 
            singleGene_transcript_donorSNP_ranges[transcript] = pd.DataFrame(cur_transcript_donor_ranges, columns=['ensembl_exon_id','donor_range_start', 'donor_range_end'])

        # finally, store these for the gene overall
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

    
    # =========================================================
    # STEP 5: keep only SNPs that fall in universal/shared regions
    # =========================================================
    gene_to_chrom = pc.drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")["chromosome_name"].to_dict()
    genes=[]
    chroms=[]
    num_common_vars_in_donor_regions=[]
    final_snp_info={}
    summary_dist_info={}
    for gene, df in final_gene_objects.items():
        cur_df = final_gene_objects[gene]
        df_non_empty = cur_df[cur_df['snps'].apply(lambda x:len(x) > 0)]
        if len(df_non_empty)>0:
            # flatten the lists in the data frame
            flat_df=explode_parallel_lists(df_non_empty)
            # get the unique snp values across all of the lists
            unique_snps = set(flat_df['snp'].values)
            # for each snp, check if it's in the shared regions
            cur_shared_regions = set(universal_donor_snp_regions_byGene.get(gene, [])) # get current gene's shared regions
            if len(cur_shared_regions)==0:
                # note that there are no shared regions
                final_snp_info[gene]=['No shared donor regions']
                genes.append(gene)
                chroms.append(gene_to_chrom[gene])
                num_common_vars_in_donor_regions.append(0)
                continue
            cur_gene_snp_info=[]
            cur_cv_info=[]
            snp_group = flat_df.groupby("snp", sort=False)
            for snp in unique_snps:
                if snp_in_regions(snp, cur_shared_regions):
                    # store information on this snp
                    cur_snp_af = snp_group.get_group(snp)['snp_af'].iloc[0]
                    cur_snp_dist = set(snp_group.get_group(snp)['snp_distance'])
                    cur_snp_be = set(snp_group.get_group(snp)['editor'])
                    # cur_snp_af=flat_df[flat_df['snp']==snp]['snp_af'].values[0]
                    # cur_snp_dist=set(flat_df[flat_df['snp']==snp]['snp_distance'])
                    # cur_snp_be=set(flat_df[flat_df['snp']==snp]['editor'])
                    if len(cur_snp_be) == 2:
                        cur_snp_be = "both ABE and CBE"
                    else:
                        cur_snp_be = list(cur_snp_be)[0]
                    # calculate average snp distance from exon boundary, since it can be different in different transcripts
                    avg_snp_dist = sum(cur_snp_dist)/len(cur_snp_dist)
                    cur_gene_snp_info.append([snp,cur_snp_af,cur_snp_dist,avg_snp_dist, cur_snp_be])
                    cur_cv_info.append([snp,cur_snp_af])
            summary_dist_info[gene]=cur_gene_snp_info
            final_snp_info[gene]=cur_cv_info
            genes.append(gene)
            chroms.append(gene_to_chrom[gene])
            num_common_vars_in_donor_regions.append(len(cur_cv_info))
        else:
            final_snp_info[gene]=[]
            genes.append(gene)
            chroms.append(gene_to_chrom[gene])
            num_common_vars_in_donor_regions.append(0)

    # convert into proper formats for saving
    be_summary_df = pd.DataFrame([
            [gene, snp, af, dist, avg_dist, editor]
            for gene, entries in summary_dist_info.items()
            for snp, af, dist,avg_dist, editor in entries],
        columns=["gene", "snp", "af", "distance_to_exon_boundary", "avg_dist", "editor"])
    donor_cv_df = pd.DataFrame({
        'gene':genes,
        'num_common_vars_in_donor_regions':num_common_vars_in_donor_regions,
        'chrom':chroms
    })
    # save the information here
    donor_cv_df.to_csv(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_summary.txt" ,sep='\t')
    donor_cv_df.to_csv(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" ,sep='\t', index=False, header=False)
    with open(output_dir + "/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl",'wb') as file:
        pickle.dump(final_snp_info,file)
    be_summary_df.to_csv(output_dir + "/ubiq_region_CommonVars/base_editor_summary.txt", sep='\t', index=False)




# --------------------------------
if __name__ == '__main__':
    main()
