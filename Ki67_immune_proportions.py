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

celltypes = ["B_cell",
    "CD4_T",
    "CD8_T",
    "Treg",
    "M1_like_macrophage",
    "M2_like_macrophage"]

pretty_labels = {
    "B_cell": "B-cells",
    "CD4_T": "CD4+ T cells",
    "CD8_T": "CD8+ T cells",
    "Treg": "Regulatory T-cells",
    "M1_like_macrophage": "M1-like macrophages",
    "M2_like_macrophage": "M2-like macrophages"}

min_cells_per_group = 20

df = ad.obs[["response_EV","cell_type_plot","Ki67_local_90_positive"]].copy()

df["response_EV"] = df["response_EV"].astype(str)
df = df[df["response_EV"].isin(order)]
df = df[df["cell_type_plot"].isin(celltypes)]

df["Ki67_local_90_positive"] = pd.to_numeric(df["Ki67_local_90_positive"],errors="coerce")

df = df.dropna(subset=["Ki67_local_90_positive"])
df["Ki67_local_90_positive"] = df["Ki67_local_90_positive"].astype(int)

summary = (df.groupby(["response_EV", "cell_type_plot"]).agg(total_cells=("Ki67_local_90_positive", "size"),ki67_positive_cells=("Ki67_local_90_positive", "sum")).reset_index())
summary["proportion_ki67_positive"] = (summary["ki67_positive_cells"] / summary["total_cells"])
summary = summary[summary["total_cells"] >= min_cells_per_group]
summary["group"] = "other"
summary.loc[summary["response_EV"].str.startswith("PR"), "group"] = "Responder"
summary.loc[summary["response_EV"].str.startswith("PD"), "group"] = "Non-responder"
summary = summary[summary["group"].isin(["Responder", "Non-responder"])]

positions = []
data = []
pvals = {}

pos = 1

for ct in celltypes:
    r = summary[
        (summary["cell_type_plot"] == ct) &
        (summary["group"] == "Responder")
    ]["proportion_ki67_positive"]

    nr = summary[
        (summary["cell_type_plot"] == ct) &
        (summary["group"] == "Non-responder")
    ]["proportion_ki67_positive"]

    data.extend([r, nr])
    positions.extend([pos, pos + 0.6])

    if len(r) > 0 and len(nr) > 0:
        pvals[ct] = mannwhitneyu(r, nr, alternative="two-sided")[1]
    else:
        pvals[ct] = np.nan

    print(f"{pretty_labels[ct]}: p = {pvals[ct]:.4f}")

    pos += 2

plt.figure(figsize=(13, 5))

bp = plt.boxplot(data,
    positions=positions,
    widths=0.45,
    patch_artist=True)

color_responder = "#4c72b0"
color_nonresponder = "#dd8452"

for i, box in enumerate(bp["boxes"]):
    if i % 2 == 0:
        box.set_facecolor(color_responder)
    else:
        box.set_facecolor(color_nonresponder)

    box.set_alpha(0.4)
    box.set_zorder(1)

for element in ["whiskers", "caps", "medians"]:
    for item in bp[element]:
        item.set_zorder(1)

for i, (x_pos, y) in enumerate(zip(positions, data)):
    if i % 2 == 0:
        c = color_responder
    else:
        c = color_nonresponder
    x = np.random.normal(x_pos, 0.04, size=len(y))
    plt.scatter(x,y,color=c,alpha=0.85,s=25,zorder=3)

for i, ct in enumerate(celltypes):
    mid = (positions[i * 2] + positions[i * 2 + 1]) / 2
    label = "p=NA" if np.isnan(pvals[ct]) else f"p={pvals[ct]:.3f}"
    plt.text(mid, 1.04, label, ha="center", fontsize=9)

group_centers = [
    (positions[i * 2] + positions[i * 2 + 1]) / 2
    for i in range(len(celltypes))]

plt.xticks(group_centers,
    [pretty_labels[ct] for ct in celltypes],
    rotation=35,
    ha="right")

plt.ylabel("Proportion of Ki67+ cells")
plt.title("Ki67+ immune cell proportions in responders vs non-responders")
plt.ylim(0, 1.1)
plt.scatter([], [], color=color_responder, label="Responder")
plt.scatter([], [], color=color_nonresponder, label="Non-responder")
plt.legend()
plt.tight_layout()
plt.show()
