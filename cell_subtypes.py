import scanpy as sc

ad = sc.read_h5ad("anndata_TissUUmaps.h5ad")

required_cols = ["immune_general",
    "tumor_general",
    "unassigned_general",
    "CD45_local_90_positive",
    "CD3e_local_90_positive",
    "CD4_local_90_positive",
    "CD8_local_90_positive",
    "FOXP3_local_90_positive",
    "CD20_local_90_positive",
    "CD68_local_90_positive",
    "CD163_local_90_positive"]

missing = [c for c in required_cols if c not in ad.obs.columns]
if missing:
    raise ValueError(f"Missing required columns in ad.obs: {missing}")

immune_mask = ad.obs["immune_general"] == 1
tumor_mask_general = ad.obs["tumor_general"] == 1

cd45  = ad.obs["CD45_local_90_positive"] == 1
cd3e  = ad.obs["CD3e_local_90_positive"] == 1
cd4   = ad.obs["CD4_local_90_positive"] == 1
cd8   = ad.obs["CD8_local_90_positive"] == 1
foxp3 = ad.obs["FOXP3_local_90_positive"] == 1
cd20  = ad.obs["CD20_local_90_positive"] == 1
cd68  = ad.obs["CD68_local_90_positive"] == 1
cd163 = ad.obs["CD163_local_90_positive"] == 1

ad.obs["Treg_rule"] = (immune_mask & cd45 & cd3e & cd4 & foxp3).astype(int)
ad.obs["CD4_T_rule"] = (immune_mask & cd45 & cd3e & cd4 & (~foxp3)).astype(int)
ad.obs["CD8_T_rule"] = (immune_mask & cd45 & cd3e & cd8).astype(int)
ad.obs["B_cell_rule"] = (immune_mask & cd45 & cd20).astype(int)
ad.obs["M2_like_macrophage_rule"] = (immune_mask & cd45 & cd68 & cd163).astype(int)
ad.obs["M1_like_macrophage_rule"] = (immune_mask & cd45 & cd68 & (~cd163)).astype(int)

rule_cols = [
    "Treg_rule",
    "CD4_T_rule",
    "CD8_T_rule",
    "B_cell_rule",
    "M1_like_macrophage_rule",
    "M2_like_macrophage_rule"
]

ad.obs["phenotype_number"] = ad.obs[rule_cols].sum(axis=1)

ad.obs["cell_type_sub"] = "non_immune"

ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 0), "cell_type_sub"] = "unassigned"

ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] > 1), "cell_type_sub"] = "ambiguous"

ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["Treg_rule"] == 1), "cell_type_sub"] = "Treg"
ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["CD4_T_rule"] == 1), "cell_type_sub"] = "CD4_T"
ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["CD8_T_rule"] == 1), "cell_type_sub"] = "CD8_T"
ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["B_cell_rule"] == 1), "cell_type_sub"] = "B_cell"
ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["M1_like_macrophage_rule"] == 1), "cell_type_sub"] = "M1_like_macrophage"
ad.obs.loc[immune_mask & (ad.obs["phenotype_number"] == 1) & (ad.obs["M2_like_macrophage_rule"] == 1), "cell_type_sub"] = "M2_like_macrophage"

ad.obs["cell_type_plot"] = "unassigned"

ad.obs.loc[tumor_mask_general, "cell_type_plot"] = "tumor"

ad.obs.loc[ad.obs["cell_type_sub"] == "Treg", "cell_type_plot"] = "Treg"
ad.obs.loc[ad.obs["cell_type_sub"] == "CD4_T", "cell_type_plot"] = "CD4_T"
ad.obs.loc[ad.obs["cell_type_sub"] == "CD8_T", "cell_type_plot"] = "CD8_T"
ad.obs.loc[ad.obs["cell_type_sub"] == "B_cell", "cell_type_plot"] = "B_cell"
ad.obs.loc[ad.obs["cell_type_sub"] == "M1_like_macrophage", "cell_type_plot"] = "M1_like_macrophage"
ad.obs.loc[ad.obs["cell_type_sub"] == "M2_like_macrophage", "cell_type_plot"] = "M2_like_macrophage"

ad.obs.loc[ad.obs["cell_type_sub"] == "ambiguous", "cell_type_plot"] = "ambiguous"

# all remaining cells stay labeled as "unassigned"

ad.obs["tumor_plot"] = (ad.obs["cell_type_plot"] == "tumor").astype(int)
ad.obs["Treg_plot"] = (ad.obs["cell_type_plot"] == "Treg").astype(int)
ad.obs["CD4_T_plot"] = (ad.obs["cell_type_plot"] == "CD4_T").astype(int)
ad.obs["CD8_T_plot"] = (ad.obs["cell_type_plot"] == "CD8_T").astype(int)
ad.obs["B_cell_plot"] = (ad.obs["cell_type_plot"] == "B_cell").astype(int)
ad.obs["M1_like_macrophage_plot"] = (ad.obs["cell_type_plot"] == "M1_like_macrophage").astype(int)
ad.obs["M2_like_macrophage_plot"] = (ad.obs["cell_type_plot"] == "M2_like_macrophage").astype(int)
ad.obs["ambiguous_plot"] = (ad.obs["cell_type_plot"] == "ambiguous").astype(int)
ad.obs["unassigned_plot"] = (ad.obs["cell_type_plot"] == "unassigned").astype(int)

ad.obs["Treg_sub"] = (ad.obs["cell_type_sub"] == "Treg").astype(int)
ad.obs["CD4_T_sub"] = (ad.obs["cell_type_sub"] == "CD4_T").astype(int)
ad.obs["CD8_T_sub"] = (ad.obs["cell_type_sub"] == "CD8_T").astype(int)
ad.obs["B_cell_sub"] = (ad.obs["cell_type_sub"] == "B_cell").astype(int)
ad.obs["M1_like_macrophage_sub"] = (ad.obs["cell_type_sub"] == "M1_like_macrophage").astype(int)
ad.obs["M2_like_macrophage_sub"] = (ad.obs["cell_type_sub"] == "M2_like_macrophage").astype(int)
ad.obs["ambiguous_sub"] = (ad.obs["cell_type_sub"] == "ambiguous").astype(int)
ad.obs["unassigned_sub"] = (ad.obs["cell_type_sub"] == "unassigned").astype(int)

print("\ncell_type_sub counts:")
print(ad.obs["cell_type_sub"].value_counts(dropna=False))

print("\ncell_type_plot counts:")
print(ad.obs["cell_type_plot"].value_counts(dropna=False))

print("\nphenotype_number counts:")
print(ad.obs["phenotype_number"].value_counts(dropna=False).sort_index())

# check that each cell has only one final assigned subtype
final_sub_cols = [
    "Treg_sub",
    "CD4_T_sub",
    "CD8_T_sub",
    "B_cell_sub",
    "M1_like_macrophage_sub",
    "M2_like_macrophage_sub",
    "ambiguous_sub",
    "unassigned_sub"
]

immune_final_sum = ad.obs.loc[immune_mask, final_sub_cols].sum(axis=1)
print("\nImmune cells with invalid subtype:")
print((immune_final_sum != 1).sum())

plot_cols = [
    "tumor_plot",
    "Treg_plot",
    "CD4_T_plot",
    "CD8_T_plot",
    "B_cell_plot",
    "M1_like_macrophage_plot",
    "M2_like_macrophage_plot",
    "ambiguous_plot",
    "unassigned_plot"
]

plot_sum = ad.obs[plot_cols].sum(axis=1)

ad.write("anndata_TissUUmaps.h5ad")
