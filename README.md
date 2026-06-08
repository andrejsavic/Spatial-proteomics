Full pipeline starting with annotation of artifacts and ending in neighborhood analysis.

1. remove_artifacts.py
Open QPTIFF in QuPath, followed by annotation of artifacts using MultiPolygon tool. The artifacts should be exported as a GeoJSON file (in this code 'newroi.geojson').
remove_artifacts.py will remove all cells whose centroid coordinates (x,y) fall within the annotated areas.

2. artifact_removal_check.py
Used to control whether the removal of artifacts using remove_artifacts.py was successful. The result is an overlay of the QPTIFF along with outlined areas where the
annotations are. The areas within the colored annotations should be white, indicating that the cells have been removed from the dataframe.

3. multiotsu.py
After running PIPEX segmentation > generate_geojson > analysis > generate_tissuumaps, the main output is anndata_TissUUmaps.h5ad. multiotsu.py allows the user to manually
decide the desired threshold for each marker, paired with visualization in TissUUmaps, allows for confident marker positivity thresholding before phenotyping cells.

4. general_celltype.py
A first-line, general phenotyping of cells. Tumor cells are in this dataset defined by being Pan-Cytokeratin positive, while immune cells by being CD45 positive. Note that
cells that are positive for both Pan-Cytokeratin AND CD45 are labeled as immune rather than tumor or ambiguous. The reason is that immune cells infilitrated in tumor areas
might also be positive for Pan-CK due to imperfect segmentation. Forcing such cells to be labeled as immune cells allows us to later study immune infiltration.

5. cell_subtypes.py
Subtyping of immune cells. Depending on markers used, the labels can be changed or added. Cells that are positive for enough markers to be labeled as multiple cell types
are labeled as ambiguous, while cells that do not meed marker positivity requirements for any cell type are labeled unassigned.

6. cell_count_plot.py
A vertical bar chart showing counts of segmented cells per patient. Note that it does not contain information about assigned/ambiguous/unassigned cells.

7. assigned_cell_proportions.py
A stacked bar chart showing proportions of assigned, ambiguous and unassigned cell proportions per patient.

8. immune_tumor_proportions.py
A stacked bar chart showing fractions of immune respectively tumor cells out of the cells assigned using general_celltype.py. Second part of the code features a
Mann-Whitney U test comparison of differences between responders and non-responders.

10. immune_subtype_proportions.py
A stacked bar chart showing proportions of immune subtypes per patient. The code includes CD8+ T-cells, CD4+ T-cells, regulatory T-cells, B-cells, M1-like macrophages
as well as M2-like macrophages but can easily be adjusted for further cell types. Contains Mann-Whitney U-test for comparing differences between responders and non-responders.

11. ADCheatmap.py
Heatmap showing proportions of tumor (Pan-Cytokeratin+) cells that are also positive for each of the three antibody-drug conjugate targets (TROP2, HER2, Nectin-4).

12. Nec4_violin.py
Violin plot showing normalized per-cell expression (signal intensity) of Nectin-4 for each patient.

13. Ki67_tumor_proportions.py
Proportions of tumor cells that are also positive for proliferation marker Ki67, along with a Mann-Whitney U test for comparing responder/non-responder differences.

14. PD1_Tcell_proportions.py
Stacked bar chart showing proportions of CD8+ and CD4+ T-cells that are positive (and negative) for the marker PD-1. Positive cells have a darker color while negative have lighter.

15. Ki67_immune_proportions.py
Boxplot comparing diffences between patients in fractions of immune subtypes that are proliferating (Ki67+).

16. distance_analysis.py
Distance analysis with a modifiable threshold. 'cell_type_from' and 'cell_type_to' as well as their markers can easily be changed depending on which two cell types you are interested in.

17. Neighborhood enrichment analysis
NE analysis looking at 10 nearest neighbors, set to plot 'tumor-tumor', 'tumor-immune', 'immune-immune' and further top 5 most statistically significant interactions.
19. 
