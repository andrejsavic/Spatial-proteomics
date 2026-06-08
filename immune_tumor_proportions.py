import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

order = ["PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"]

df = ad.obs[["response_EV", "immune_general", "tumor_general"]].copy()

df["general_type_plot"] = "unassigned"
df.loc[df["immune_general"] == 1, "general_type_plot"] = "immune"
df.loc[df["tumor_general"] == 1, "general_type_plot"] = "tumor"

df = df[df["general_type_plot"].isin(["immune", "tumor"])]
df = df[df["response_EV"].isin(order)]

count_df = (df.groupby(["response_EV", "general_type_plot"]).size().reset_index(name="count"))

count_df["proportion"] = (count_df.groupby("response_EV")["count"].transform(lambda x: x / x.sum()))

plot_df = count_df.pivot(index="response_EV",
    columns="general_type_plot",
    values="proportion").fillna(0)

plot_df = plot_df.reindex(order).fillna(0)
plot_df = plot_df.reindex(columns=["immune", "tumor"]).fillna(0)
plot_df = plot_df.rename(columns={"immune": "CD45+ immune cells",
    "tumor": "PanCK+ tumor cells"})

stats_df = plot_df.reset_index()

stats_df["response_group"] = "Other"
stats_df.loc[stats_df["response_EV"].str.startswith("PR"), "response_group"] = "Responder"
stats_df.loc[stats_df["response_EV"].str.startswith("PD"), "response_group"] = "Non-responder"

r = stats_df[stats_df["response_group"] == "Responder"]["CD45+ immune cells"].dropna()
nr = stats_df[stats_df["response_group"] == "Non-responder"]["CD45+ immune cells"].dropna()

stat, pval = mannwhitneyu(r, nr, alternative="two-sided")

print(f"CD45+ immune cells: Mann-Whitney U p-value = {pval:.4f}")

plot_df = plot_df.iloc[::-1]

ax = plot_df.plot(kind="barh",
    stacked=True,
    figsize=(12, 6),
    width=0.85)

ax.set_xlim(0, 1)
ax.set_xlabel("Proportion of cells")
ax.set_ylabel("Patients")
ax.set_title("CD45+ immune and PanCK+ tumor cell proportions per patient")
ax.legend(title="Cell category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left")
plt.tight_layout()
plt.show()

box_data = stats_df[stats_df["response_group"].isin(["Responder", "Non-responder"])].copy()
group_order = ["Responder", "Non-responder"]

palette = {
    "Responder": "#1f77b4",
    "Non-responder": "#ff7f0e"
}

plt.figure(figsize=(5, 6))

ax = sns.boxplot(data=box_data,
    x="response_group",
    y="CD45+ immune cells",
    order=group_order,
    palette=palette,
    showfliers=False,
    width=0.5)

sns.stripplot(data=box_data,
    x="response_group",
    y="CD45+ immune cells",
    order=group_order,
    palette=palette,
    alpha=0.8,
    linewidth=0.5,
    size=6,
    jitter=0.12,
    ax=ax)

y_min = box_data["CD45+ immune cells"].min()
y_max = box_data["CD45+ immune cells"].max()
y_range = y_max - y_min

if y_range == 0:
    y_range = 0.1

ax.text(0.5,
    y_max + 0.08 * y_range,
    f"p={pval:.3f}",
    ha="center")

plt.ylabel("Proportion of CD45+ immune cells")
plt.xlabel("Response group")
plt.title("CD45+ immune cell proportions per patient")
plt.ylim(max(0, y_min - 0.05 * y_range),
    min(1, y_max + 0.20 * y_range))
plt.tight_layout()
plt.show()
