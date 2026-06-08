import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

celltypes_to_plot = ["Treg",
    "CD4_T",
    "CD8_T",
    "B_cell",
    "M1_like_macrophage",
    "M2_like_macrophage"]

order = ["PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"]

df = ad.obs[["response_EV", "cell_type_plot"]].copy()

df = df[df["cell_type_plot"].isin(celltypes_to_plot)]
df = df[df["response_EV"].isin(order)]

count_df = (df.groupby(["response_EV", "cell_type_plot"]).size().reset_index(name="count"))
count_df["proportion"] = count_df.groupby("response_EV")["count"].transform(lambda x: x / x.sum())

plot_df = count_df.pivot(index="response_EV",
    columns="cell_type_plot",
    values="proportion").fillna(0)

plot_df = plot_df.reindex(order).fillna(0)
plot_df = plot_df.reindex(columns=celltypes_to_plot).fillna(0)
plot_df = plot_df.rename(columns={
    "Treg": "Regulatory T-cells",
    "CD4_T": "CD4+ T cells",
    "CD8_T": "CD8+ T cells",
    "B_cell": "B-cells",
    "M1_like_macrophage": "M1-like macrophages",
    "M2_like_macrophage": "M2-like macrophages"})
plot_df = plot_df.iloc[::-1]

ax = plot_df.plot(kind="barh",
    stacked=True,
    figsize=(12, 6),
    width=0.85)

ax.set_xlim(0, 1)
ax.set_xlabel("Proportion of cells")
ax.set_ylabel("Patients")
ax.set_title("Immune cell type proportions per patient")
ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
plt.show()

stats_df = plot_df.reset_index()
stats_df["group"] = "other"
stats_df.loc[stats_df["response_EV"].str.startswith("PR"), "group"] = "Responder"
stats_df.loc[stats_df["response_EV"].str.startswith("PD"), "group"] = "Non-responder"
celltype_cols = [
    "Regulatory T-cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "B-cells",
    "M1-like macrophages",
    "M2-like macrophages"]

plot_long = stats_df.melt(
    id_vars=["response_EV", "group"],
    value_vars=celltype_cols,
    var_name="cell_type",
    value_name="proportion")

plot_long = plot_long[plot_long["group"].isin(["Responder", "Non-responder"])]

positions = []
data = []
pvals = {}

pos = 1
for celltype in celltype_cols:
    r = plot_long[
        (plot_long["cell_type"] == celltype) &
        (plot_long["group"] == "Responder")
    ]["proportion"]

    nr = plot_long[
        (plot_long["cell_type"] == celltype) &
        (plot_long["group"] == "Non-responder")
    ]["proportion"]
    data.extend([r, nr])
    positions.extend([pos, pos + 0.6])
    pvals[celltype] = mannwhitneyu(r, nr, alternative="two-sided")[1]
    print(f"{celltype}: Mann-Whitney U p-value = {pvals[celltype]:.4f}")
    pos += 2

plt.figure(figsize=(13, 5))

plt.boxplot(data, positions=positions, widths=0.45)

for x_pos, y in zip(positions, data):
    x = np.random.normal(x_pos, 0.04, size=len(y))
    plt.scatter(x, y)

for i, celltype in enumerate(celltype_cols):
    mid = (positions[i * 2] + positions[i * 2 + 1]) / 2
    plt.text(
        mid,
        1.04,
        f"p={pvals[celltype]:.3f}",
        ha="center",
        fontsize=9)

group_centers = [
    (positions[i * 2] + positions[i * 2 + 1]) / 2
    for i in range(len(celltype_cols))]

plt.xticks(group_centers, celltype_cols, rotation=35, ha="right")
plt.ylabel("Proportion of cells")
plt.title("Immune cell type proportions in responders vs non-responders")
plt.ylim(0, 1.1)

plt.tight_layout()
plt.show()
