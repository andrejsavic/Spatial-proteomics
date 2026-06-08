from skimage.filters import threshold_otsu, threshold_multiotsu
import seaborn as sns
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad("anndata_TissUUmaps.h5ad")

df = adata.to_df()

markers = [
    "Nectin4_local_90",
    "PD-1_local_90",
    "CD163_local_90",
    "FOXP3_local_90",
    "HER2_local_90",
    "Ki67_local_90",
    "CD8_local_90",
    "CD20_local_90",
    "TROP2_local_90",
    "PDL1_local_90",
    "Pan-Cytokeratin_local_90",
    "CD68_local_90",
    "CD3e_local_90",
    "CD4_local_90",
    "CD45_local_90"
]

def get_threshold(array, title="", method="otsu", plot_histogram=True):
    array = np.asarray(array).flatten()
    array = array[np.isfinite(array)]

    if method == "otsu":
        threshold = threshold_otsu(array)
    elif method == "otsu1":
        thresholds = threshold_multiotsu(array, classes=3)
        threshold = thresholds[0]
    elif method == "otsu2":
        thresholds = threshold_multiotsu(array, classes=3)
        threshold = thresholds[1]
    elif method == "otsu3":
        thresholds = threshold_multiotsu(array, classes=4)
        threshold = thresholds[2]
    elif method == "otsu4":
        thresholds = threshold_multiotsu(array, classes=5)
        threshold = thresholds[3]
    elif method == "otsu5":
        thresholds = threshold_multiotsu(array, classes=6)
        threshold = thresholds[4]
    else:
        raise ValueError(f"Unknown thresholding method: {method}")

#    if plot_histogram:
#        p = sns.displot(np.arcsinh(array), kind="kde", height=3)
#        p.fig.set_dpi(100)
#        ymax = p.ax.get_ylim()[1]
#        plt.plot(
#            [np.arcsinh(threshold), np.arcsinh(threshold)],
#            [0, ymax],
#            color="red"
#        )
#        plt.title(f"{title} histogram with {method}")
#        plt.xlabel("arcsinh(expression)")
#        plt.ylabel("Density")
#        plt.show()
#
    return threshold

marker_methods = {
    "Nectin4_local_90": "otsu2",
    "PD-1_local_90": "otsu3",
    "CD163_local_90": "otsu4",
    "FOXP3_local_90": "otsu4",
    "HER2_local_90": "otsu2",
    "Ki67_local_90": "otsu3",
    "CD8_local_90": "otsu4",
    "CD20_local_90": "otsu5",
    "TROP2_local_90": "otsu",
    "PDL1_local_90": "otsu4",
    "Pan-Cytokeratin_local_90": "otsu2",
    "CD68_local_90": "otsu3",
    "CD3e_local_90": "otsu3",
    "CD4_local_90": "otsu4",
    "CD45_local_90": "otsu2",
}

threshold_summary = {}

for marker in markers:
    print("\n" + "=" * 70)
    method = marker_methods.get(marker, "otsu")
    print(f"Marker: {marker}")
    print(f"Method: {method}")

    array = df[marker].to_numpy()
    threshold = get_threshold(array, title=marker, method=method, plot_histogram=True)

    positive_col = f"{marker}_positive"
    adata.obs[positive_col] = (df[marker] > threshold).astype(int)

    percent_positive = adata.obs[positive_col].mean() * 100
    number_positive = int(adata.obs[positive_col].sum())

    print(f"Threshold = {threshold}")
    print(f"Percent positive = {percent_positive:.2f}%")
    print(f"Number positive = {number_positive}")

    threshold_summary[marker] = {
        "method": method,
        "threshold": float(threshold),
        "percent_positive": float(percent_positive),
        "number_positive": number_positive
    }

print("\n" + "#" * 70)
print("SUMMARY")
for marker, info in threshold_summary.items():
    print(
        f"{marker}: "
        f"method={info['method']}, "
        f"threshold={info['threshold']:.6f}, "
        f"percent_positive={info['percent_positive']:.2f}%, "
        f"number_positive={info['number_positive']}"
    )

adata.write_h5ad("anndata_TissUUmaps.h5ad")
