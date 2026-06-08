import json
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import shape
from shapely.ops import unary_union

df = pd.read_csv("cell_data_artifacts_removed.csv")

with open("newroi.geojson", "r") as f:
    gj = json.load(f)

geoms = [shape(feat["geometry"]) for feat in gj.get("features", [])]
roi = unary_union(geoms)

def plot_roi(g):
    if g.geom_type == "Polygon":
        x, y = g.exterior.xy
        plt.plot(x, y, linewidth=2)
    elif g.geom_type == "MultiPolygon":
        for p in g.geoms:
            x, y = p.exterior.xy
            plt.plot(x, y, linewidth=2)

plt.figure(figsize=(7, 7))
plt.scatter(df["x"], df["y"], s=1)
plot_roi(roi)
plt.gca().set_aspect("equal", adjustable="box")
plt.gca().invert_yaxis()
plt.title("Filtered cells + ROI")
plt.xlabel("x"); plt.ylabel("y")
plt.show()
