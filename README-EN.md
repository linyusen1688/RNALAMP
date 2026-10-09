# RNALAMP: Spatial RNA Colocalization Community Analysis

## 1. Project Purpose
xxx


## 2. Prerequisites

### 2.1 Python Dependencies

The example code uses or depends on the following Python packages:
- `spatialdata`
- `spatialdata-plot`
- `anndata`
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `networkx`
- `igraph`
- `leidenalg`
- `statsmodels`
- `shapely`
- `scikit-image`
- `tqdm`

Run the code in a Python environment that can access the project data and the `rnalamp` source code, and make sure all the packages listed above are installed.

## Calculating CLScore
```python
python cal_clscore.py --input /path_to_point_file --output /path_to_save_file --tmp-dir /path_to_save_temp_data
```
The program outputs a table with the following columns:

| Column | Field | Description / Processing |
|---|---|---|
| 1 | `observe_gene` | Observed gene |
| 2 | `reference_gene` | Reference gene |
| 3 | `p_val` | Colocalization test p-value |
| 4 | `min_p_value` | Minimum p-value |
| 5 | `effect_area` | Effect area |
| 6 | `last_intersection` | Stored as `None` if the value is `None`; otherwise converted to a float |
| 7 | `max_enrichment` | Maximum enrichment |
| 8 | `max_z_score` | Maximum Z-score |
| 9 | `max_enrichment_radius` | Radius at which maximum enrichment occurs |
| 10 | `max_z_radius` | Radius at which the maximum Z-score occurs |

CLScore can also be calculated for selected cell types by setting the `--select-cell` argument. Note that annotations for these cell types must be present in the input file.
```python
python cal_clscore.py --input /path_to_point_file --output /path_to_save_dir/gene_colocation.csv --tmp-dir /path_to_save_temp_data --select-cell Gra Basal-1 Basal-2 Spi-1 Spi-2 
```
## Loading the Data

Import the required packages
```python
import spatialdata_plot  # Register the SpatialData plotting extension
import spatialdata as sd
from rnalamp import rnalamp
```
```python
zarr_path = "zhenfeng.AFL.spRNAseq.P1-NC.zarr"
sdata = sd.read_zarr(zarr_path)
sdata
```

```python
SpatialData object, with associated Zarr store: /mnt/dfc_data3/project/linyusen/database/07_spatialdata/output-XETG00546__0076988__TMA__20251120__071742/annonate3/zhenfeng.AFL.spRNAseq.P1-NC.zarr
├── Images
│     └── 'morphology_focus': DataTree[cyx] (4, 11797, 12000), (4, 5898, 6000), (4, 2949, 3000), (4, 1475, 1500), (4, 737, 750)
├── Labels
│     ├── 'cell_labels': DataTree[yx] (11797, 12000), (5898, 6000), (2949, 3000), (1475, 1500), (737, 750)
│     └── 'nucleus_labels': DataTree[yx] (11797, 12000), (5898, 6000), (2949, 3000), (1475, 1500), (737, 750)
├── Points
│     └── 'transcripts': DataFrame with shape: (<Delayed>, 13) (3D points)
├── Shapes
│     ├── 'cell_boundaries': GeoDataFrame shape: (2603, 1) (2D shapes)
│     ├── 'cell_circles': GeoDataFrame shape: (2603, 2) (2D shapes)
│     └── 'nucleus_boundaries': GeoDataFrame shape: (2478, 1) (2D shapes)
└── Tables
      └── 'table': AnnData (2322, 4491)
with coordinate systems:
    ▸ 'global', with elements:
        morphology_focus (Images), cell_labels (Labels), nucleus_labels (Labels), transcripts (Points), cell_boundaries (Shapes), cell_circles (Shapes), nucleus_boundaries (Shapes)
```
This dataset contains one large microscopy image (`morphology_focus`), as well as label images representing the segmentation masks of cells and nuclei (`cell_labels` and `nucleus_labels`, respectively).


Next, convert the `SpatialData` object into a format usable by `rnalamp` for downstream analysis. The `rnalamp` class inherits from `SpatialData`, so `SpatialData` methods can also be called on `rnalamp` objects.
```python
rdata = rnalamp(images=sdata.images,
                labels=sdata.labels,
                points=sdata.points,
                shapes=sdata.shapes,
                tables=sdata.tables)
```

Load the gene-pair colocalization results calculated by `cal_clscore.py`.
```python
clscore_file='/path_to_save_dir/gene_colocation.csv'
rdata.read_clscore(input_file=clscore_file) #读取计算好的结果
```


## Colocalization Analysis
In this step, we construct a gene colocalization network using the CLScore values and significance-test results from the gene colocalization analysis. We then cluster the network with the Leiden community detection algorithm, grouping genes with significant colocalization relationships into relatively cohesive functional communities. These communities provide the basis for subsequent community-quality assessment and regional enrichment analysis.
```python
rdata.build_gene_colocation_communities(pval_threshold=0.01,resolution=1.5,random_seed=42)  # Set parameters for community detection
```

After detecting communities in the gene colocalization network, we perform quality control to assess the validity of the identified community structures from two main perspectives:

**Local Modularity Contribution:** A permutation test compares the local modularity contribution of each observed community with the distribution obtained by randomly shuffling node-to-community assignments. This assesses whether the community's network modular structure is statistically significant.
**Conductance:** A permutation test compares the conductance of each observed community with the corresponding distribution for randomized communities, assessing the relationship between within-community and between-community connections. In general, lower conductance indicates that connections within a community are denser relative to connections crossing its boundary.

The parameter `n_permutations=500` specifies 500 random permutations for each test. The program outputs the statistics for each community and the number of communities that pass the significance threshold of `p < 0.05` set in the code.
Community QC results can help determine whether the current community partition is reasonable. If many communities are very small or lack cohesive structure, consider lowering `resolution` in light of the QC results. Conversely, increasing it may reveal a finer-grained community structure. Parameter tuning should take network topology, the distribution of community sizes, and biological context into account rather than relying on a single metric.
```python
rdata.community_qc(n_permutations=500)  # Quality-check the detected communities; assess whether an overly high resolution creates communities that are too small
```

Finally, we use a set of organelle-localization marker genes compiled from APEX-seq data to perform regional enrichment analysis on each detected gene colocalization community. This evaluates the overall organelle-localization tendencies of different communities and helps explore the spatial functional characteristics of the gene colocalization network.
```python
gene_loc_file = "/mnt/dfc_data3/project/linyusen/database/17_RNA_colocation/01_marker/data/gene_loc.csv" # Organelle-localization marker gene set
rdata.community_enrich_region(gene_loc_file=gene_loc_file) # Calculate the cellular regions enriched in each community
```

## Saving Results
After community detection, export each gene colocalization network community separately and save it in a text format readable by Cytoscape.

One `community_<community_ID>.txt` file is generated for each community. Each file contains a tab-separated edge list with the following fields:
- `source`: Source gene of the edge
- `target`: Target gene of the edge
- `p_val`: Colocalization p-value
- `effect_area`: Effect area
These files can be used for subsequent network inspection or imported into Cytoscape.
```python
output_dir = '/path_to_savedir/'
rdata.save_communities(output_dir=output_dir) # Save communities in a format that can be imported into Cytoscape
```
