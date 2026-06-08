import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

ad = sc.read_h5ad("merged_anndata_with_response.h5ad")

cell_type_from = "CD4_T"
cell_type_to = "CD8_T"
marker_from = None
marker_to = None
distance_um = 20
pixel_size_um = 0.5077

df = ad.obs.copy()
df = df[df["cell_type_plot"].isin([cell_type_from, cell_type_to])]

res = []

for p in df["response_EV"].unique():
    d = df[df["response_EV"] == p]

    from_df = d[d["cell_type_plot"] == cell_type_from]
    to_df = d[d["cell_type_plot"] == cell_type_to]

    if marker_from:
        from_df = from_df[from_df[marker_from] == 1]
    if marker_to:
        to_df = to_df[to_df[marker_to] == 1]

    if len(from_df) == 0 or len(to_df) == 0:
        continue

    tree = cKDTree(to_df[["x","y"]].values)
    dist = tree.query(from_df[["x","y"]].values, k=1)[0] * pixel_size_um

    frac = (dist <= distance_um).mean()

    res.append({"response_EV": p, "close": frac, "far": 1-frac})

res = pd.DataFrame(res)

order = ["PR1","PR2","PR3","PR4","PR5","PR6","PR7","PR8","PR9","SD","PD1","PD2","PD3","PD4","PD5","N/K1","N/K2","N/K3"]
res["response_EV"] = pd.Categorical(res["response_EV"], categories=order, ordered=True)
res = res.sort_values("response_EV")

plt.figure(figsize=(10,7))
plt.barh(res["response_EV"].astype(str), res["close"], label="≤20 µm")
plt.barh(res["response_EV"].astype(str), res["far"], left=res["close"], label=">20 µm")
plt.xlim(0,1)
plt.xlabel(f"Proportion of PanCK+ T-cells")
plt.title(f"PanCK+ cells proximity to M2-like macrophages ({distance_um} µm)")
plt.legend()
plt.tight_layout()
plt.show()

res["group"] = "other"
res.loc[res["response_EV"].str.startswith("PR"), "group"] = "Responder"
res.loc[res["response_EV"].str.startswith("PD"), "group"] = "Non-responder"

r = res[res["group"]=="Responder"]["close"]
nr = res[res["group"]=="Non-responder"]["close"]

pval = mannwhitneyu(r, nr)[1]
print("p-value:", pval)

plt.figure(figsize=(5,5))
plt.boxplot([r, nr], tick_labels=["Responder","Non-responder"])

for i, y in enumerate([r, nr], 1):
    x = np.random.normal(i, 0.04, size=len(y))
    plt.scatter(x, y)

plt.ylabel(f"PanCK+ cells within {distance_um} µm from M2-like macrophages")
plt.title(f"Responder vs Non-responder\np = {pval:.3e}")
plt.ylim(0,1)
plt.tight_layout()
plt.show()
