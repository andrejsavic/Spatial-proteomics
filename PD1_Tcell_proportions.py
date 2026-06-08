import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

order = [
    "PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"
]

df = ad.obs[[
    "response_EV",
    "cell_type_plot",
    "PD-1_local_90_positive"
]].copy()

df["response_EV"] = df["response_EV"].astype(str)
df = df[df["response_EV"].isin(order)]
df = df[df["cell_type_plot"].isin(["CD4_T", "CD8_T"])]

df["category"] = ""

df.loc[
    (df["cell_type_plot"] == "CD8_T") & (df["PD-1_local_90_positive"] == 1),
    "category"
] = "CD8+ T cells PD-1+"

df.loc[
    (df["cell_type_plot"] == "CD4_T") & (df["PD-1_local_90_positive"] == 1),
    "category"
] = "CD4+ T cells PD-1+"

df.loc[
    (df["cell_type_plot"] == "CD8_T") & (df["PD-1_local_90_positive"] == 0),
    "category"
] = "CD8+ T cells PD-1-"

df.loc[
    (df["cell_type_plot"] == "CD4_T") & (df["PD-1_local_90_positive"] == 0),
    "category"
] = "CD4+ T cells PD-1-"

plot_order = [
    "CD8+ T cells PD-1+",
    "CD4+ T cells PD-1+",
    "CD8+ T cells PD-1-",
    "CD4+ T cells PD-1-"
]

colors = [
    "#c44e52",  # medium red (CD8 PD1+)
    "#4c72b0",  # medium blue (CD4 PD1+)
    "#e89a9a",  # lighter red (CD8 PD1-)
    "#9eb6e6"   # lighter blue (CD4 PD1-)
]

count_df = (
    df.groupby(["response_EV", "category"])
      .size()
      .reset_index(name="count")
)

count_df["proportion"] = count_df.groupby("response_EV")["count"].transform(lambda x: x / x.sum())

plot_df = count_df.pivot(
    index="response_EV",
    columns="category",
    values="proportion"
).fillna(0)

plot_df = plot_df.reindex(order).fillna(0)
plot_df = plot_df.reindex(columns=plot_order).fillna(0)

ax = plot_df.plot(
    kind="barh",
    stacked=True,
    figsize=(12, 7),
    width=0.85,
    color=colors
)

ax.set_xlim(0, 1)
ax.set_xlabel("Proportion of total CD4+ and CD8+ T cells")
ax.set_ylabel("Patients")
ax.set_title("PD-1 positivity among CD4+ and CD8+ T cells per patient")
ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
plt.show()

# ---- stats on CD8 PD1+ ----
stats_df = plot_df.reset_index()

stats_df["group"] = "other"
stats_df.loc[stats_df["response_EV"].str.startswith("PR"), "group"] = "Responder"
stats_df.loc[stats_df["response_EV"].str.startswith("PD"), "group"] = "Non-responder"

r = stats_df[stats_df["group"] == "Responder"]["CD8+ T cells PD-1+"].dropna()
nr = stats_df[stats_df["group"] == "Non-responder"]["CD8+ T cells PD-1+"].dropna()

pval = mannwhitneyu(r, nr, alternative="two-sided")[1]
print(f"CD8+ T cells PD-1+: Mann-Whitney U p-value = {pval:.4f}")

plt.figure(figsize=(5, 5))
plt.boxplot([r, nr], tick_labels=["Responder", "Non-responder"])

for i, y in enumerate([r, nr], 1):
    x = np.random.normal(i, 0.04, size=len(y))
    plt.scatter(x, y)

plt.ylabel("Proportion of total CD4+ and CD8+ T cells")
plt.title(f"PD-1+ CD8+ T-cell proportions per patient\np = {pval:.3e}")
plt.ylim(0, 1)

plt.tight_layout()
plt.show()
