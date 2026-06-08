import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

expr = ad.to_df()
meta = ad.obs.copy()

df = meta.join(expr)

df = df[df["cell_type_plot"] == "tumor"].copy()
df = df[df["Nectin4_local_90"].notna()].copy()

order = ["PR1", "PR2", "PR3", "PR4", "PR5", "PR6", "PR7", "PR8", "PR9",
    "SD",
    "PD1", "PD2", "PD3", "PD4",
    "N/K1", "N/K2", "N/K3"]

df["response_EV"] = df["response_EV"].astype(str)
df = df[df["response_EV"].isin(order)]

patient_order = order

data_to_plot = [
    df.loc[df["response_EV"] == pid, "Nectin4_local_90"].values
    for pid in patient_order]

fig, ax = plt.subplots(figsize=(12, 6))

parts = ax.violinplot(
    data_to_plot,
    showmeans=False,
    showmedians=True,
    showextrema=False)

for i, pid in enumerate(patient_order, start=1):
    y = df.loc[df["response_EV"] == pid, "Nectin4_local_90"].values
    x = np.random.normal(i, 0.06, size=len(y))
    ax.scatter(x, y, s=5, alpha=0.2)

ax.set_xticks(range(1, len(patient_order) + 1))
ax.set_xticklabels(patient_order, rotation=45, ha="right")
ax.set_xlabel("Patients")
ax.set_ylabel("Nectin4 expression (local_90)")
ax.set_title("Nectin4 expression in tumor cells per patient")

plt.tight_layout()
plt.show()
