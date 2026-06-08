import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

order = ["PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"]

df = ad.obs.loc[
    ad.obs["cell_type_plot"] == "tumor",
    ["response_EV",
     "TROP2_local_90_positive",
     "HER2_local_90_positive",
     "Nectin4_local_90_positive"]].copy()

df = df[df["response_EV"].isin(order)]

prop_df = (
    df.groupby("response_EV")[[
        "TROP2_local_90_positive",
        "HER2_local_90_positive",
        "Nectin4_local_90_positive"]].mean())

prop_df = prop_df.rename(columns={
    "TROP2_local_90_positive": "TROP2",
    "HER2_local_90_positive": "HER2",
    "Nectin4_local_90_positive": "Nectin4"})

prop_df = prop_df[["Nectin4", "TROP2", "HER2"]]
prop_df = prop_df.reindex(order)

plt.figure(figsize=(6, 10))
plt.imshow(prop_df, aspect="auto")
plt.colorbar(label="Proportion positive")

plt.xticks(range(len(prop_df.columns)), prop_df.columns)
plt.yticks(range(len(prop_df.index)), prop_df.index)

plt.xlabel("Marker")
plt.ylabel("Patients")
plt.title("Marker positivity in PanCK+ cells per patient")

plt.tight_layout()
plt.show()

stats_df = prop_df.reset_index()

stats_df["group"] = "other"
stats_df.loc[stats_df["response_EV"].str.startswith("PR"), "group"] = "Responder"
stats_df.loc[stats_df["response_EV"].str.startswith("PD"), "group"] = "Non-responder"

r = stats_df[stats_df["group"] == "Responder"]["Nectin4"].dropna()
nr = stats_df[stats_df["group"] == "Non-responder"]["Nectin4"].dropna()

pval = mannwhitneyu(r, nr, alternative="two-sided")[1]

print(f"Nectin4+ PanCK+ cells: Mann-Whitney U p-value = {pval:.4f}")

plt.figure(figsize=(5, 5))
plt.boxplot([r, nr], tick_labels=["Responder", "Non-responder"])

for i, y in enumerate([r, nr], 1):
    x = np.random.normal(i, 0.04, size=len(y))
    plt.scatter(x, y)

plt.ylabel("Proportion of cells")
plt.title(f"Nectin4+ PanCK+ cell proportions per patient\np = {pval:.3e}")
plt.ylim(0, 1)

plt.tight_layout()
plt.show()
