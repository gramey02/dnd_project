# imports
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory to save files to.')
    parser.add_argument('--exon_file',type=str,required=True,help="Filename for list of chromosomes and genes of interest.")
    parser.add_argument('--nearest_gene_file',type=str,required=True,help="File that contains information on what the coordinates of genes throughout the genome are.")
    parser.add_argument('--upstream_excision_dist',type=int,required=True,help="distance upstream of the gene start to look for variants.")
    parser.add_argument('--downstream_excision_dist', type=int,required=True,help="distance downstream of the gene end to look for variants.")
    parser.add_argument('--filter_out_nearby_genes', type=int, required=True, help="Boolean--filter out nearby protein coding genes from excision window.")
    args = parser.parse_args()
    return args

def restrict_window(goi_start, goi_end, window_start, window_end, other_genes):
    """
    Restrict the GOI ± flank window to avoid overlaps with other genes.
    The restricted window may exclude parts of the GOI if necessary.

    Parameters
    ----------
    goi_start, goi_end : int
        Start and end coordinates of the gene of interest (GOI).
    flank : int
        Bases to extend on each side.
    other_genes : list of (start, end)
        Coordinates of other genes on the same chromosome.

    Returns
    -------
    tuple or None
        (restricted_start, restricted_end) if any overlap with GOI remains,
        otherwise None.
    """

    # sort and merge any overlapping other_gene intervals
    other_genes = sorted(other_genes, key=lambda x: x[0])
    merged = []
    for s, e in other_genes:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)

    # Start with full window, subtract all other genes
    intervals = [(window_start, window_end)]
    for s, e in merged:
        new_intervals = []
        for a, b in intervals:
            if e < a or s > b:
                # no overlap
                new_intervals.append((a, b))
            else:
                # trim overlapping part
                if s > a:
                    new_intervals.append((a, s - 1))
                if e < b:
                    new_intervals.append((e + 1, b))
        intervals = new_intervals

    # Keep only subranges that still overlap the GOI
    goi_intervals = []
    for a, b in intervals:
        if b >= goi_start and a <= goi_end:
            goi_intervals.append((a, b))

    if not goi_intervals:
        return None  # no overlap with GOI left

    # Return the largest subrange that overlaps GOI
    restricted = max(goi_intervals, key=lambda x: x[1] - x[0])
    return restricted

def overlaps(a, b):
    """Return True if two ranges (start, end) overlap at all."""
    return not (a[1] < b[0] or b[1] < a[0])

def encompasses(a,b):
    """Return True if range a is fully encompassed by range b."""
    return (a[0]>=b[0] and a[1]<=b[1])

def encompasses_key_exon(excision_range,gene,exon_df):
    # filter exon df to get the gene's exon information
    # filter exon df to get the gene's exon information
    cur_gene_exons=exon_df[exon_df.hgnc_symbol==gene]
    # determine if, in at least one of the transcripts, the range encompasses at least one exon
    result=False
    for transcript in cur_gene_exons.ensembl_transcript_id.unique():
        transcript_df=cur_gene_exons[cur_gene_exons.ensembl_transcript_id==transcript]
        num_exons=max(transcript_df['rank'].values)

        # get ranges of exons in current transcript
        exon_ranges=[]
        for idx,row in transcript_df.iterrows():
            exon_ranges.append((row.exon_chrom_start,row.exon_chrom_end))

        # for each exon range, 
        for exon_range in exon_ranges:
            if encompasses(exon_range,excision_range):
                result=True
                return result
    return result

def main():
    # load input args
    args=parse_args()
    output_dir=args.output_dir
    exon_df=pd.read_csv(args.exon_file, index_col=0)
    genebody_file=args.nearest_gene_file
    upstream_excision_dist=args.upstream_excision_dist
    downstream_excision_dist=args.downstream_excision_dist
    filter_out_nearby_genes=args.filter_out_nearby_genes
    #----------------
    exon_filt=exon_df[['hgnc_symbol','chromosome_name','start_position','end_position','strand','ensembl_gene_id']].drop_duplicates()
    # set default excision coords
    exon_filt_minus=exon_filt[exon_filt.strand==-1]
    exon_filt_plus=exon_filt[exon_filt.strand==1]
    # add new columns for excision start and end to plus strand genes
    exon_filt_plus = exon_filt_plus.assign(excision_start=exon_filt_plus['start_position']-upstream_excision_dist)
    exon_filt_plus = exon_filt_plus.assign(excision_end=exon_filt_plus['end_position']+downstream_excision_dist)

    # add new columns for start and end to minus strand genes
    exon_filt_minus = exon_filt_minus.assign(excision_start=exon_filt_minus['start_position']-downstream_excision_dist)
    exon_filt_minus = exon_filt_minus.assign(excision_end=exon_filt_minus['end_position'] + upstream_excision_dist)

    # remerge data frames
    excision_df1=pd.concat([exon_filt_minus,exon_filt_plus])

    # ----------------
    if filter_out_nearby_genes==1:
        # load file that contains coordinates for various gene bodies
        cols = ["chrom","source","feature","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(
            genebody_file,   # works with .gz too
            sep="\t",
            comment="#",
            names=cols,
            dtype={
                "seqname": "string",
                "source": "string",
                "feature": "category",
                "start": "int64",
                "end": "int64",
                "score": "string",    # '.' or a number → keep as string
                "strand": "category",
                "frame": "string"     # '.', '0', '1', '2'
            },
            keep_default_na=False,    # keep '.' as literal, not NaN
            engine="c",               # default; fast
        )
        # parse out the attributes in the file
        def parse_gtf_attributes(attr: str) -> dict:
            out = {}
            # 1) Trim whitespace and the trailing ';', then split into fields
            for field in attr.strip().rstrip(";").split(";"):
                field = field.strip()
                if not field:
                    continue
                # 2) GTF form: key "value"
                if " " in field:
                    k, v = field.split(" ", 1)        # split once: left is key, right is `"value"`
                    out[k] = v.strip().strip('"')     # remove surrounding quotes
                # 3) Also accept GFF3 form: key=value
                elif "=" in field:
                    k, v = field.split("=", 1)
                    out[k] = v
            return out
        attr_df = pd.DataFrame(df["attribute"].map(parse_gtf_attributes).tolist())
        df = pd.concat([df.drop(columns="attribute"), attr_df], axis=1)

        # conduct some filtering on the file
        ref_chroms=['chr'+str(x) for x in list(range(1,23))] # valid ref chromosomes
        ref_chroms=ref_chroms+['chrX']
        df_filt=df[df['chrom'].isin(ref_chroms)]
        # removing version number from ensembl gene id
        df_filt[['gene_id_base', 'gene_id_ver']] = (
        df_filt['gene_id'].astype('string')
        .str.split('.', n=1, expand=True, regex=False)   # use regex=False or escape the dot
        )

        # exclude genes of interest from gene body file
        df_pc=df_filt[df_filt.gene_type=='protein_coding']
        df_external_genes=df_pc[~df_pc['gene_name'].isin(excision_df1.hgnc_symbol)]
        df_external_genes=df_external_genes[~df_external_genes['gene_id_base'].isin(excision_df1.ensembl_gene_id)]
        df_external_genes=df_external_genes[df_external_genes['feature']=='gene']
        # create a dictionary that stores genebodies on each chromosome (reduces filtering steps/computational time later)
        filt_df_dict={}
        for chrom in df_external_genes.chrom.unique():
            df_chrom=df_external_genes[df_external_genes.chrom==chrom]
            filt_df_dict[chrom]=df_chrom

        # now, for each gene and its excision coordinates, see if there is a protein-coding gene overlapping, and restrict its coordinates accordingly
        min_excision_window=[]
        max_excision_window=[]
        good_excision_candidate=[]
        for idx,row in excision_df1.iterrows():
            # filter the genebody df to the current chromosome
            excision_start=row.excision_start
            excision_end=row.excision_end
            gene_start=row.start_position
            gene_end=row.end_position
            cur_genebody_df = filt_df_dict['chr'+str(row.chromosome_name)]

            # filter gene bodies to those that overlap with the excision window
            gene_ranges = []
            for idx_range,row_range in cur_genebody_df.iterrows():
                gene_ranges.append((row_range.start,row_range.end))
            window = (excision_start, excision_end)
            overlapping_genes = [r for r in gene_ranges if overlaps(r, window)] # overlapping gene ranges

            # get the restricted range of excision based on overlapping genes
            result=restrict_window(gene_start,gene_end,excision_start,excision_end,overlapping_genes)

            if result is not None:
                # we can also check if the window encompasses the start or end exon
                #print(row.hgnc_symbol)
                good_for_excision=encompasses_key_exon(result,row.hgnc_symbol,exon_df)
                # save this gene's bounds
                min_excision_window.append(result[0])
                max_excision_window.append(result[1])
                good_excision_candidate.append(good_for_excision)
            else:
                min_excision_window.append(None)
                max_excision_window.append(None)
                good_excision_candidate.append(False)

        # combine into df
        excision_df2=pd.DataFrame({
            'hgnc_symbol':excision_df1.hgnc_symbol,
            'chrom':excision_df1.chromosome_name,
            'excision_start':excision_df1.excision_start,
            'excision_end':excision_df1.excision_end,
            'new_excision_start':min_excision_window,
            'new_excision_end':max_excision_window,
            'gene_start':excision_df1.start_position,
            'gene_end':excision_df1.end_position,
            'strand':excision_df1.strand,
            'ensembl_gene_id':excision_df1.ensembl_gene_id,
            'good_excision_candidate':good_excision_candidate
        })
        excision_df2['excision_window_size']=excision_df2['new_excision_end']-excision_df2['new_excision_start']

    else:
        excision_df2=pd.DataFrame({
            'hgnc_symbol':excision_df1.hgnc_symbol,
            'chrom':excision_df1.chromosome_name,
            'excision_start':excision_df1.excision_start,
            'excision_end':excision_df1.excision_end,
            'new_excision_start':excision_df1.excision_start,
            'new_excision_end':excision_df1.excision_end,
            'gene_start':excision_df1.start_position,
            'gene_end':excision_df1.end_position,
            'strand':excision_df1.strand,
            'ensembl_gene_id':excision_df1.ensembl_gene_id,
            'good_excision_candidate':1
        })
        excision_df2['excision_window_size']=excision_df2['new_excision_end']-excision_df2['new_excision_start']

    
    # save file
    excision_df2.to_csv(output_dir + "/gene_filtered_excision_coords.txt",sep='\t')
    excision_df2.to_csv(output_dir + "/gene_filtered_excision_coords_noIDX.txt",sep='\t',header=False,index=False)





# ------------------
if __name__=='__main__':
    main()