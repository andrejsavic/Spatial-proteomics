import scanpy as sc

ad = sc.read_h5ad("anndata_TissUUmaps.h5ad")

cols_to_remove = ["general_type",
    "immune_general",
    "tumor_general",
    "unassigned_general"]

cols_existing = [c for c in cols_to_remove if c in ad.obs.columns]
ad.obs = ad.obs.drop(columns=cols_existing)

ad.obs["immune_general"] = 0
ad.obs["tumor_general"] = 0
ad.obs["unassigned_general"] = 0

immune_mask = ad.obs["CD45_local_90_positive"] == 1
ad.obs.loc[immune_mask, "immune_general"] = 1

tumor_mask = ((ad.obs["Pan-Cytokeratin_local_90_positive"] == 1) & (ad.obs["immune_general"] == 0))

ad.obs.loc[tumor_mask, "tumor_general"] = 1

unassigned_mask = ((ad.obs["immune_general"] == 0) & (ad.obs["tumor_general"] == 0))

ad.obs.loc[unassigned_mask, "unassigned_general"] = 1

group_sum = (ad.obs["immune_general"] + ad.obs["tumor_general"] + ad.obs["unassigned_general"])

print("Immune:", int(ad.obs["immune_general"].sum()))
print("Tumor:", int(ad.obs["tumor_general"].sum()))
print("Unassigned:", int(ad.obs["unassigned_general"].sum()))
print("Invalid assignments:", int((group_sum != 1).sum()))

ad.write("anndata_TissUUmaps_general_phenotypes.h5ad")
