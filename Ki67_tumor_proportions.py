import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

ad = sc.read_h5ad("merged_anndata_TissUUmaps.h5ad")

df = ad.obs.loc[
    ad.obs["cell_type_plot"] == "tumor",
    ["patientid", "Ki67_local_90_positive"]
].copy()

df["Ki67_status"] = df["Ki67_local_90_positive"].map({
    1: "Ki67+",
    0: "Ki67-"
})

count_df = (
    df.groupby(["patientid", "Ki67_status"])
      .size()
      .reset_index(name="count")
)

count_df["proportion"] = count_df.groupby("patientid")["count"].transform(lambda x: x / x.sum())
plot_df = count_df.pivot(index="patientid", columns="Ki67_status", values="proportion").fillna(0)
plot_df = plot_df.reindex(columns=["Ki67+", "Ki67-"])

try:
    plot_df.index = plot_df.index.astype(int)
    plot_df = plot_df.sort_index()
    plot_df.index = plot_df.index.astype(str)
except:
    plot_df = plot_df.sort_index()

plot_df = plot_df.iloc[::-1]

ax = plot_df.plot(
    kind="barh",
    stacked=True,
    figsize=(8, 10),
    width=0.85
)

ax.set_xlabel("Proportion of PanCK+ cells")
ax.set_ylabel("Patient ID")
ax.set_title("Ki67 positivity among PanCK+ cells")
ax.set_xlim(0, 1)
ax.legend(title="Ki67 status", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()
