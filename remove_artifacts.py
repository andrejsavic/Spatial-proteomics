import json
import pandas as pd
from shapely.geometry import shape, Point
from shapely.ops import unary_union

df = pd.read_csv("cell_data.csv")

with open("newroi.geojson", "r") as f:
    gj = json.load(f)

geoms = []
for feat in gj["features"]:
    g = shape(feat["geometry"])
    if g.geom_type in ["LineString", "MultiLineString"]:
        g = g.buffer(100)
    geoms.append(g)

roi_union = unary_union(geoms)

inside = df.apply(lambda r: roi_union.contains(Point(r["x"], r["y"])),
    axis=1)

print("Cells excluded:", inside.sum())

df_filtered = df.loc[~inside].copy()
df_filtered.to_csv("cell_data_artifacts_removed.csv", index=False)
