#!/usr/bin/env python3
"""Format the master D&D summary dataframe for downstream use."""

import argparse
import pickle
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


DATE = datetime.now().strftime("%Y-%m-%d")

SCRIPT_DIR = Path(__file__).resolve().parent
DND_PROJECT_ROOT = SCRIPT_DIR.parent.parent

TARGET_COLUMN_ORDER = [
    "hgnc_symbol",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "MOI",
    "CLASSIFICATION",
    "HI Score",
    "%HI",
    "pLI",
    "LOEUF",
    "ensg",
    "chrom",
    "obs_lof",
    "exp_lof",
    "prior_mean",
    "post_mean",
    "post_lower_95",
    "post_upper_95",
    "indel_targetable_post_NMD_assessment",
    "num_vars_inducing_NMD",
    "crisproff_targetable",
    "num_crisproff_vars",
    "base_editable",
    "num_base_editable_vars",
    "excision_targetable",
    "num_excision_vars",
    "num_excision_pairs",
    "hets_across_strats_prePAM",
    "prop_hets_across_strats_prePAM",
    "num_indel_hets_prePAM",
    "indel_hets_prop_prePAM",
    "num_crisproff_hets_prePAM",
    "crisproff_hets_prop_prePAM",
    "num_base_edit_hets_prePAM",
    "base_edit_hets_prop_prePAM",
    "num_excision_hets_prePAM",
    "excision_hets_prop_prePAM",
    "indel_pam_targetable",
    "crisproff_pam_targetable",
    "base_edit_pam_targetable",
    "excision_pam_targetable",
    "hets_across_strats",
    "prop_hets_across_strats",
    "num_indel_hets",
    "indel_hets_prop",
    "num_crisproff_hets",
    "crisproff_hets_prop",
    "num_base_edit_hets",
    "base_edit_hets_prop",
    "num_excision_hets",
    "excision_hets_prop",
    "indels_num_haps_four_guides_alg2",
    "indels_num_hets_four_guides_alg2",
    "indels_num_haps_all_guides_alg2",
    "indels_num_hets_all_guides_alg2",
    "indels_max_guides",
    "CRISPRoff_num_haps_four_guides_alg2",
    "CRISPRoff_num_hets_four_guides_alg2",
    "CRISPRoff_num_haps_all_guides_alg2",
    "CRISPRoff_num_hets_all_guides_alg2",
    "CRISPRoff_max_guides",
    "combined_base_edit_num_haps_four_guides_alg2",
    "combined_base_edit_num_hets_four_guides_alg2",
    "combined_base_edit_num_haps_all_guides_alg2",
    "combined_base_edit_num_hets_all_guides_alg2",
    "combined_base_edit_max_guides",
    "excision_num_haps_four_guides_alg2",
    "excision_num_hets_four_guides_alg2",
    "excision_num_haps_all_guides_alg2",
    "excision_num_hets_all_guides_alg2",
    "excision_max_guides",
]

RENAME_MAP = {
    "post_mean": "s_het",
    "post_lower_95": "s_het_lower_95",
    "post_upper_95": "s_het_upper_95",
    "indel_targetable_post_NMD_assessment": "exon_disr_targetable",
    "num_vars_inducing_NMD": "num_exon_disr_vars",
    "crisproff_targetable": "epi_sil_targetable",
    "num_crisproff_vars": "num_epi_sil_vars",
    "base_editable": "ss_disr_targetable",
    "num_base_editable_vars": "num_ss_disr_vars",
    "indel_pam_targetable": "exon_disr_spcas9_targetable",
    "crisproff_pam_targetable": "epi_sil_spcas9_targetable",
    "base_edit_pam_targetable": "ss_disr_spcas9_targetable",
    "excision_pam_targetable": "excision_spcas9_targetable",
    "hets_across_strats": "num_hets_all_strats",
    "prop_hets_across_strats": "prop_1KG_hets_all_strats",
    "num_indel_hets": "num_hets_exon_disr",
    "indel_hets_prop": "prop_1KG_hets_exon_disr",
    "num_crisproff_hets": "num_hets_epi_sil",
    "crisproff_hets_prop": "prop_1KG_hets_epi_sil",
    "num_base_edit_hets": "num_hets_ss_disr",
    "base_edit_hets_prop": "prop_1KG_hets_ss_disr",
    "num_excision_hets": "num_hets_excision",
    "excision_hets_prop": "prop_1KG_hets_excision",
    "hets_across_strats_prePAM": "num_hets_preCas_all_strats",
    "prop_hets_across_strats_prePAM": "prop_hets_preCas_all_strats",
    "num_indel_hets_prePAM": "num_hets_preCas_exon_disr",
    "indel_hets_prop_prePAM": "prop_hets_preCas_exon_disr",
    "num_crisproff_hets_prePAM": "num_hets_preCas_epi_sil",
    "crisproff_hets_prop_prePAM": "prop_hets_preCas_epi_sil",
    "num_base_edit_hets_prePAM": "num_hets_preCas_ss_disr",
    "base_edit_hets_prop_prePAM": "prop_hets_preCas_ss_disr",
    "num_excision_hets_prePAM": "num_hets_preCas_excision",
    "excision_hets_prop_prePAM": "prop_hets_preCas_excision",
    "indels_num_haps_four_guides_alg2": "num_haplos_targeted_4g_exon_disr",
    "indels_num_hets_four_guides_alg2": "num_hets_targeted_4g_exon_disr",
    "indels_num_haps_all_guides_alg2": "num_haplos_targeted_ag_exon_disr",
    "indels_num_hets_all_guides_alg2": "num_hets_targeted_ag_exon_disr",
    "indels_max_guides": "max_num_guides_exon_disr",
    "CRISPRoff_num_haps_four_guides_alg2": "num_haplos_targeted_4g_epi_sil",
    "CRISPRoff_num_hets_four_guides_alg2": "num_hets_targeted_4g_epi_sil",
    "CRISPRoff_num_haps_all_guides_alg2": "num_haplos_targeted_ag_epi_sil",
    "CRISPRoff_num_hets_all_guides_alg2": "num_hets_targeted_ag_epi_sil",
    "CRISPRoff_max_guides": "max_num_guides_epi_sil",
    "combined_base_edit_num_haps_four_guides_alg2": "num_haplos_targeted_4g_ss_disr",
    "combined_base_edit_num_hets_four_guides_alg2": "num_hets_targeted_4g_ss_disr",
    "combined_base_edit_num_haps_all_guides_alg2": "num_haplos_targeted_ag_ss_disr",
    "combined_base_edit_num_hets_all_guides_alg2": "num_hets_targeted_ag_ss_disr",
    "combined_base_edit_max_guides": "max_num_guides_ss_disr",
    "excision_num_haps_four_guides_alg2": "num_haplos_targeted_4g_excision",
    "excision_num_hets_four_guides_alg2": "num_hets_targeted_4g_excision",
    "excision_num_haps_all_guides_alg2": "num_haplos_targeted_ag_excision",
    "excision_num_hets_all_guides_alg2": "num_hets_targeted_ag_excision",
    "excision_max_guides": "max_num_guides_excision",
}

TARGETABLE_COLS = [
    "exon_disr_targetable",
    "epi_sil_targetable",
    "ss_disr_targetable",
    "excision_targetable",
]

SPCAS9_TARGETABLE_COLS = [
    "exon_disr_spcas9_targetable",
    "epi_sil_spcas9_targetable",
    "ss_disr_spcas9_targetable",
    "excision_spcas9_targetable",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format the master D&D summary dataframe for downstream analyses."
    )
    parser.add_argument(
        "--input-csv",
        required=True,
        help="CSV produced by Creating_Master_DnD_DataFrame.py.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="Optional output path. Defaults to <input_stem>_formatted_<DATE>.csv.",
    )
    parser.add_argument(
        "--hpo-pkl",
        required=True,
        help="Pickle mapping HGNC symbols to HPO term summaries.",
    )
    parser.add_argument(
        "--dominant-mutations-csv",
        required=True,
        help="ClinVar/OMIM dominant mutation count CSV.",
    )
    return parser.parse_args()


def validate_columns(df: pd.DataFrame, required_cols: list[str], label: str) -> None:
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing {label} columns: {', '.join(missing)}")


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    validate_columns(df, TARGET_COLUMN_ORDER, "input")
    return df[TARGET_COLUMN_ORDER].copy()


def add_gene_group(df: pd.DataFrame) -> pd.DataFrame:
    df.insert(1, "Gene_Group", np.where(df["s_het"] < 0.1, "Group 2", "Group 1"))
    return df


def add_hpo_terms(df: pd.DataFrame, hpo_pkl: str) -> pd.DataFrame:
    with open(hpo_pkl, "rb") as fp:
        gene_to_hpo = pickle.load(fp)
    df.insert(7, "HPO_terms", df["hgnc_symbol"].map(gene_to_hpo))
    return df


def add_dominant_mutation_counts(df: pd.DataFrame, dominant_mutations_csv: str) -> pd.DataFrame:
    dom_clean = pd.read_csv(dominant_mutations_csv)
    if "gene" in dom_clean.columns:
        dom_clean = dom_clean.rename(columns={"gene": "hgnc_symbol"})

    required_cols = ["hgnc_symbol", "dominant_mutation_count"]
    validate_columns(dom_clean, required_cols, "dominant mutation")
    df = df.merge(dom_clean, on="hgnc_symbol", how="left")
    dominant_col = df.pop("dominant_mutation_count")
    df.insert(7, "dominant_mutation_count", dominant_col)
    return df


def fill_targetable_columns(df: pd.DataFrame) -> pd.DataFrame:
    df[TARGETABLE_COLS] = df[TARGETABLE_COLS].fillna(0)
    df[SPCAS9_TARGETABLE_COLS] = df[SPCAS9_TARGETABLE_COLS].fillna(0)
    return df


def build_output_path(input_csv: str, output_csv: str | None) -> Path:
    if output_csv:
        return Path(output_csv)

    input_path = Path(input_csv)
    return input_path.with_name(f"{input_path.stem}_formatted_{DATE}.csv")


def main() -> None:
    args = parse_args()

    data = pd.read_csv(args.input_csv)
    data = reorder_columns(data)
    data = data.rename(columns=RENAME_MAP)
    data = data[data["CLASSIFICATION"].isin(["Strong", "Definitive", "Moderate"])].copy()
    data = add_gene_group(data)
    data = add_hpo_terms(data, args.hpo_pkl)
    data = add_dominant_mutation_counts(data, args.dominant_mutations_csv)

    if "counter" in data.columns:
        data = data.drop(columns=["counter"])

    data = fill_targetable_columns(data)

    output_path = build_output_path(args.input_csv, args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
