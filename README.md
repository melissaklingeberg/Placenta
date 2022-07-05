Immunohistochemical analysis of cell-type specific genes during placental development

This script was created in the context of a 5-month student project at The Human Protein Atlas, with the aim to analyse placenta cell-type specific genes by tissue microarray-based immunohistochemistry to identify trophoblast specific marker. This code can be used as a basis to analyse annotation data from other experiments, identify candidates of interest and present the results graphically.

Description:
In section two, a figure of the TMA was generated. Afterwards, quality control was performed for the tissue cores (section 3) and antibodies (section 4). Tissue cores, which did not include at least one cell type of interest in at least 25 % of all slides were excluded. For antibody QC, the values for the staining intensity were categorized into undistinct (negative, weak) and distinct (moderate, strong). Antibodies were excluded, which did not show distinct staining in at least two tissue cores. Both approaches were presented graphically. Finally, filtering methods were applied to screen for potential marker proteins for trophoblasts in general and trophoblast sub types (section 5). Trophoblast marker and their subtype specificity were visualized in a Venn Diagram and an Upset Plot, their expression over placental development in a connected scatterplot (section 6). 

Get started:
Annotation Raw data parameters incorporate:
-	Patient ID 
-	Timepoint of pregnancy
-	Antibody ID
-	different cell types 
-	staining intensity (0 - negative, 1 - weak, 2 - moderate, 3 - strong)

The file „Genes“ gives information about the protein-specificity of the antibody. 
The file “Patient Information” contains more detailed information about patients incl their position on the TMA.
The file “tropho_spec_manual” includes a list of the final marker proteins (manual reviewing was needed due to lack in annotation data).
