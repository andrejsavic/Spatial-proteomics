import scanpy as sc
import matplotlib.pyplot as plt

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

assigned_types = ["tumor",
    "Treg",
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

df["assignment_status"] = "unassigned"
df.loc[df["cell_type_plot"].isin(assigned_types), "assignment_status"] = "assigned"
df.loc[df["cell_type_plot"] == "ambiguous", "assignment_status"] = "ambiguous"

count_df = (df.groupby(["response_EV", "assignment_status"]).size().reset_index(name="count"))
count_df["proportion"] = count_df.groupby("response_EV")["count"].transform(lambda x: x / x.sum())

plot_df = count_df.pivot(
    index="response_EV",
    columns="assignment_status",
    values="proportion"
).fillna(0)

plot_df = plot_df.reindex(order)
plot_df = plot_df.reindex(columns=["assigned", "ambiguous", "unassigned"])

plot_df = plot_df.rename(columns={
    "assigned": "Assigned",
    "ambiguous": "Ambiguous",
    "unassigned": "Unassigned"
})

plot_df = plot_df.iloc[::-1]

ax = plot_df.plot(kind="barh",
    stacked=True,
    figsize=(12, 6),
    width=0.85)

ax.set_xlim(0, 1)
ax.set_xlabel("Proportion of cells")
ax.set_ylabel("EV response")
ax.set_title("Assigned, unassigned and ambiguous cells by EV response")
ax.legend(title="Status", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()
