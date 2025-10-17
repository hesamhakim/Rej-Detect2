import pandas as pd


def preprocess_variant_tables(
    recipient_file: str,
    donor_file: str,
    mixed_file: str,
    mean_depth: int = 800,
    donor_homoplasmy_threshold: float = 0.95,
    recipient_low_threshold: float = 0.05,
):
    df_R = pd.read_csv(recipient_file, sep="\t")
    df_D = pd.read_csv(donor_file, sep="\t")
    df_RDM = pd.read_csv(mixed_file, sep="\t")

    required_cols = {"Pos", "Ref", "Variant", "VariantLevel"}
    for df in (df_R, df_D, df_RDM):
        if not required_cols.issubset(df.columns):
            missing = required_cols - set(df.columns)
            raise ValueError(f"Variant table missing required columns: {sorted(missing)}")

    for df in (df_R, df_D, df_RDM):
        df["VariantLevel"] = pd.to_numeric(df["VariantLevel"], errors="coerce").fillna(0.0)

    common_variants_df = pd.merge(df_D, df_R, on="Pos", suffixes=("_D", "_R"), how="inner")

    df_D_homoplasmic = df_D[df_D["VariantLevel"] >= donor_homoplasmy_threshold]

    merged_df = pd.merge(df_D_homoplasmic, df_R, on="Pos", suffixes=("_D", "_R"), how="left")
    merged_df["VariantLevel_R"] = merged_df["VariantLevel_R"].fillna(0.0)

    informative_positions = merged_df[merged_df["VariantLevel_R"] < recipient_low_threshold]["Pos"].unique()

    df_RDM_informative = df_RDM[df_RDM["Pos"].isin(informative_positions)].copy()

    missing_positions = set(informative_positions) - set(df_RDM_informative["Pos"])
    if missing_positions:
        df_missing = pd.DataFrame(
            {
                "Pos": list(missing_positions),
                "Ref": ["N"] * len(missing_positions),
                "Variant": ["N"] * len(missing_positions),
                "VariantLevel": [0.0] * len(missing_positions),
            }
        )
        df_RDM_informative = pd.concat([df_RDM_informative, df_missing], ignore_index=True)

    if "Coverage" in df_RDM_informative.columns:
        df_RDM_informative["n_reads"] = pd.to_numeric(df_RDM_informative["Coverage"], errors="coerce").fillna(mean_depth).astype(int)
    else:
        df_RDM_informative["n_reads"] = int(mean_depth)

    df_RDM_informative["k_reads"] = (df_RDM_informative["VariantLevel"] * df_RDM_informative["n_reads"]).round().astype(int)

    informative_df = df_RDM_informative[["Pos", "k_reads", "n_reads"]].copy().sort_values("Pos").reset_index(drop=True)

    return informative_df, common_variants_df
