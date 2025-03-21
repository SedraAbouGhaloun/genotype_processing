import argparse
import pandas as pd
import numpy as np
import os

# argparse:
file_path_parser = argparse.ArgumentParser(description=" ")
# arguments
file_path_parser.add_argument("--in_pgs_file", help="path to unzipped pgs file.")
file_path_parser.add_argument(
    "--in_pvar_extended_ids", help="path to unzipped pgs file."
)
file_path_parser.add_argument("--nvar", help="The size of the subset.", nargs="+")
file_path_parser.add_argument("--out", help="The output path.")


args = file_path_parser.parse_args()

snpefflike_file = args.in_pgs_file + ".snpefflike"

if not os.path.isfile(snpefflike_file):
    pgs = pd.read_csv(args.in_pgs_file, sep="\t", comment="#")
    pvar_df = pd.read_csv(args.in_pvar_extended_ids, sep="\t")
    if not "ID" in pvar_df.columns or "ID_x" in pvar_df.columns:
        pvar_df = pd.read_csv(args.in_pvar_extended_ids, sep=",")
    if not "ID" in pvar_df.columns:
        pvar_df = pvar_df.rename(columns={"ID_x": "ID"})
    if not "ID" in pvar_df.columns or "ID_x" in pvar_df.columns:
        raise ValueError(
            "Could not find id column in pvar file with either comma- or tab-separation."
        )

    chr_name_col = "chr_name"
    chr_pos_col = "chr_position"
    ref_allele_col = "reference_allele"
    effect_allele_col = "effect_allele"

    if not chr_name_col in pgs.columns:
        chr_name_col = "hm_chr"
    if not chr_pos_col in pgs.columns:
        chr_pos_col = "hm_pos"
    if not ref_allele_col in pgs.columns:
        if "other_allele" in pgs.columns:
            ref_allele_col = "other_allele"
        else:
            ref_allele_col = "hm_inferOtherAllele"
    if not effect_allele_col in pgs.columns:
        effect_allele_col = "hm_effectAllele"

    pgs["tmp_id"] = pgs.apply(
        lambda x: f"{x[chr_name_col]}:{x[chr_pos_col]}_{x[ref_allele_col]}_{x[effect_allele_col]}",
        axis=1,
    )

    pgs.index = pgs["tmp_id"]

    pvar_df.index = pvar_df["id"]
    pgs = pgs.join(pvar_df, rsuffix="_pvar")
    pvar_df.index = pvar_df["flipped_id"]
    pgs = pgs.join(pvar_df, rsuffix="flipped_")
    pgs.loc[~pd.notnull(pgs["ID"]), "ID"] = pgs.loc[
        ~pd.notnull(pgs["ID"]), "IDflipped_"
    ]
    pgs.index = pgs["ID"]

    print(
        f"Number of matched variants: {pd.notnull(pgs.index).sum()} (out of {len(pgs.index)}) : {int(100 * pd.notnull(pgs.index).sum() / len(pgs.index))}"
    )
    pgs = pgs[pd.notnull(pgs.index)]

    # Making it snpefflike
    snpefflike = pd.DataFrame()
    snpefflike["#CHROM"] = pgs[chr_name_col]
    snpefflike["POS"] = pgs[chr_pos_col]
    snpefflike["ID"] = pgs["ID"]
    try:
        snpefflike["REF"] = pgs["reference_allele"]
    except KeyError:
        try:
            snpefflike["REF"] = pgs["other_allele"]
        except KeyError:
            snpefflike["REF"] = pgs["hm_inferOtherAllele"]

    snpefflike["ALT"] = pgs["effect_allele"]
    snpefflike["BETA"] = pgs["effect_weight"]

    snpefflike.to_csv(snpefflike_file, sep="\t", index=False)
else:
    print(
        f"Read snpefflike file from {snpefflike_file} instead of putting it together again."
    )
    snpefflike = pd.read_csv(snpefflike_file, sep="\t")
snpefflike["ABS_BETA"] = np.abs(snpefflike["BETA"])
snpefflike = snpefflike.sort_values("ABS_BETA", ascending=False)
snpefflike = snpefflike.loc[~snpefflike["ID"].duplicated()]

betas = snpefflike["BETA"]

for suffix in args.nvar:
    n = int(suffix[:-1])
    if suffix[-1] == "k":
        n *= 1000
    elif suffix[-1] == "M":
        n *= 1000000
    print(f"Using {n} variants.")

    try:
        t = np.abs(betas[n])
        print(f"Threshold is abs(BETA) >= {t}")
    except IndexError:
        t = np.max(betas)
        print(
            f"Including all variants, as there are only {len(betas)} instead of the {n} requested ones."
        )

    subset_path = f"{args.out}/subsets_pgs/variant_subset.{suffix}.tsv"
    print(
        f'Found {(np.abs(snpefflike["BETA"].values) >= t).sum()} values over the threshold.'
    )
    subset = snpefflike.iloc[:n, :]
    print(f"Subset has shape: {subset.shape}")
    subset.to_csv(subset_path, sep="\t", index=False)

print(f"Done, wrote PRS with rsids to {subset_path}")
