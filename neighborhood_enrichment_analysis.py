import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

adata = sc.read_h5ad("merged_anndata_with_response.h5ad")

adata.obs["response_EV"] = adata.obs["response_EV"].astype(str).str.strip()
adata = adata[adata.obs["response_EV"].str.startswith(("PR", "PD"))].copy()
adata.obs["response_group"] = np.where(adata.obs["response_EV"].str.startswith("PR"), "Responder", "Non-responder")
adata = adata[adata.obs["cell_type_plot"].notna() & ~adata.obs["cell_type_plot"].isin(["ambiguous", "unassigned"])].copy()
adata.obs["cell_type_plot"] = adata.obs["cell_type_plot"].astype("category")
adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()

rows = []

for patient in adata.obs["patientid"].unique():
    ad = adata[adata.obs["patientid"] == patient].copy()

    ad.obs["cell_type_plot"] = ad.obs["cell_type_plot"].astype("category")

    sq.gr.spatial_neighbors(ad, coord_type="generic", n_neighs=10)
    np.random.seed(42)
    sq.gr.nhood_enrichment(ad, cluster_key="cell_type_plot", seed=42)

    cats = list(ad.obs["cell_type_plot"].cat.categories)

    z = pd.DataFrame(
        ad.uns["cell_type_plot_nhood_enrichment"]["zscore"],
        index=cats,
        columns=cats)

    c = pd.DataFrame(
        ad.uns["cell_type_plot_nhood_enrichment"]["count"],
        index=cats,
        columns=cats)

    z_long = (
        z.reset_index()
        .melt(id_vars="index", var_name="cell_type_2", value_name="zscore")
        .rename(columns={"index": "cell_type_1"}))

    z_long["count"] = (
        c.reset_index()
        .melt(id_vars="index", var_name="cell_type_2", value_name="count")["count"])

    z_long["patientid"] = patient
    z_long["response_EV"] = ad.obs["response_EV"].iloc[0]
    z_long["response_group"] = ad.obs["response_group"].iloc[0]

    rows.append(z_long)

results = pd.concat(rows, ignore_index=True)

results["interaction"] = results.apply(
    lambda r: "__".join(sorted([r["cell_type_1"], r["cell_type_2"]])),
    axis=1)

pairs = (
    results
    .groupby(["patientid", "response_EV", "response_group", "interaction"], as_index=False)
    .agg(zscore=("zscore", "mean"), count=("count", "mean")))

immune_celltypes = [
    "B_cell",
    "CD4_T",
    "CD8_T",
    "Treg",
    "M1_like_macrophage",
    "M2_like_macrophage"]

selected = ["tumor__tumor"] + ["__".join(sorted(["tumor", ct])) for ct in immune_celltypes]
plot_data = pairs[pairs["interaction"].isin(selected)].copy()
stats = []

for interaction in selected:
    d = pairs[pairs["interaction"] == interaction]
    r = d[d["response_group"] == "Responder"]["zscore"].dropna()
    nr = d[d["response_group"] == "Non-responder"]["zscore"].dropna()
    stat, p = mannwhitneyu(r, nr, alternative="two-sided")

    stats.append({
        "interaction": interaction,
        "p_value": p,
        "median_responder": r.median(),
        "median_non_responder": nr.median()})

stats = pd.DataFrame(stats)
print(stats)

pvals = dict(zip(stats["interaction"], stats["p_value"]))

pretty = {
    "tumor__tumor": "Tumor–tumor",
    "__".join(sorted(["tumor", "B_cell"])): "Tumor–B cell",
    "__".join(sorted(["tumor", "CD4_T"])): "Tumor–CD4+ T cell",
    "__".join(sorted(["tumor", "CD8_T"])): "Tumor–CD8+ T cell",
    "__".join(sorted(["tumor", "Treg"])): "Tumor–Treg",
    "__".join(sorted(["tumor", "M1_like_macrophage"])): "Tumor–M1-like macrophage",
    "__".join(sorted(["tumor", "M2_like_macrophage"])): "Tumor–M2-like macrophage"}

plot_data["interaction_label"] = plot_data["interaction"].map(pretty)
order = [pretty[x] for x in selected]

plt.figure(figsize=(12, 6))

ax = sns.boxplot(
    data=plot_data,
    x="interaction_label",
    y="zscore",
    hue="response_group",
    order=order,
    showfliers=False)

sns.stripplot(
    data=plot_data,
    x="interaction_label",
    y="zscore",
    hue="response_group",
    order=order,
    dodge=True,
    alpha=0.8,
    linewidth=0.5,
    ax=ax)

ax.axhline(0, color="black", linestyle="--", linewidth=1)

y_min = plot_data["zscore"].min()
y_max = plot_data["zscore"].max()
y_range = y_max - y_min

for i, interaction in enumerate(selected):
    y = plot_data[plot_data["interaction"] == interaction]["zscore"].max()
    ax.text(i, y + 0.06 * y_range, f"p={pvals[interaction]:.3f}", ha="center")

plt.xticks(rotation=45, ha="right")
plt.ylabel("Neighborhood enrichment z-score")
plt.xlabel("Cell type interaction")
plt.title("Neighborhood enrichment of tumor cells")
plt.ylim(y_min - 0.05 * y_range, y_max + 0.18 * y_range)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="Response group")

plt.tight_layout()
plt.show()
