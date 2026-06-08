import scanpy as sc
import matplotlib.pyplot as plt

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")
df = ad.obs.copy()

order = ["PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"]

count_df = (df.groupby("response_EV")
      .size()
      .reset_index(name="cell_count"))

count_df = count_df[count_df["response_EV"].isin(order)]

count_df["response_EV"] = pd.Categorical(
    count_df["response_EV"],
    categories=order,
    ordered=True)

count_df = count_df.sort_values("response_EV")
count_df["cell_count_millions"] = count_df["cell_count"] / 1000000

gap_after = ["PR9", "SD", "PD4"]

x = []
pos = 0
for sample in count_df["response_EV"]:
    x.append(pos)
    pos += 1
    if sample in gap_after:
        pos += 0.7

plt.figure(figsize=(11, 5))
plt.bar(x,count_df["cell_count_millions"])
plt.xticks(x, count_df["response_EV"], rotation=45)
plt.xlabel("Patients")
plt.ylabel("Number of segmented cells (millions)")
plt.title("Cell counts per patient")
plt.ticklabel_format(style="plain", axis="y")
plt.ylim(0, 2.5)
plt.tight_layout()
plt.show()
