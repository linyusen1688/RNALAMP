# RNALAMP：空间 RNA 共定位社区分析

## 1. 项目用途
xxx


## 2. 运行前准备

### 2.1 Python 依赖

示例代码使用或依赖以下 Python 包：
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

请在能够读取项目数据和 `rnalamp` 源码的 Python 环境中运行，并确认上述包已安装。

## 计算clscore
```python
python cal_clscore.py --input /path_to_point_file --output /path_to_save_file --tmp-dir /path_to_save_temp_data
```
该程序输出一个列表，列表每一列含义如下

| 列索引 | 字段含义 | 处理方式 |
|---|---|---|
| 1 | `observe_gene` | 观察基因 |
| 2 | `reference_gene` | 参考基因 |
| 3 | `p_val` | 共定位检验 p 值 |
| 4 | `min_p_value` | 最小 p 值 |
| 5 | `effect_area` | 效应面积 |
| 6 | `last_intersection` | 若为 `None` 则存为 `None`，否则转为浮点数 |
| 7 | `max_enrichment` | 最大富集值 |
| 8 | `max_z_score` | 最大 Z-score |
| 9 | `max_enrichment_radius` | 最大富集对应半径 |
| 10 | `max_z_radius` | 最大 Z-score 对应半径 |

我们也可以仅针对指定细胞类型计算 CLScore。通过设置参数--select-cell即可实现，需要注意的是这些细胞的注释需要存在于输入文件中
```python
python cal_clscore.py --input /path_to_point_file --output /path_to_save_dir/gene_colocation.csv --tmp-dir /path_to_save_temp_data --select-cell Gra Basal-1 Basal-2 Spi-1 Spi-2 
```
## 加载数据

导入所需的包
```python
import spatialdata_plot  # 注册 SpatialData 绘图扩展
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
这个数据集包含1张大型显微镜图像（morphology_focus），此外，它还包含标记图像，即细胞或细胞核的分割掩膜（分别为cell_labels和nucleus_labels）


接着我们需要将spatialdata转换为rnalamp可用的格式，以便于我们后续的分析，rnalamp继承自spatialdata，换而言之，spatialdata的函数均可再rnalamp中调用
```python
rdata = rnalamp(images=sdata.images,
                labels=sdata.labels,
                points=sdata.points,
                shapes=sdata.shapes,
                tables=sdata.tables)
```

导入由cal_clscore.py计算的基因对共定位结果
```python
clscore_file='/path_to_save_dir/gene_colocation.csv'
rdata.read_clscore(input_file=clscore_file) #读取计算好的结果
```


## 共定位分析
在此步骤中，我们将基于基因共定位分析得到的 clscore 及其显著性检验结果构建基因共定位网络，并采用 Leiden 社区检测算法对网络进行聚类，将具有显著共定位关系的基因划分为多个相对紧密的功能社区，为后续的社区质量评估与区域富集分析提供基础。
```python
rdata.build_gene_colocation_communities(pval_threshold=0.01,resolution=1.5,random_seed=42)  #设置参数划分社区
```

完成基因共定位网络的社区划分后，我们将对识别出的社区进行质量控制，主要从以下两个方面评估社区结构的合理性：

局部模块度贡献（Local Modularity Contribution）：通过置换检验，将真实社区的局部模块度贡献与随机打乱节点社区归属后得到的分布进行比较，评估真实社区的网络模块结构是否具有统计学上的显著性。
传导性（Conductance）：通过置换检验，将真实社区的传导性与随机社区的对应分布进行比较，评估社区内部连接与跨社区连接的相对关系。通常，传导性越低，说明社区内部连接相对于跨社区连接越紧密。

参数 n_permutations=500 表示每项检验执行 500 次随机置换。程序将输出各社区的统计结果，以及在代码设定的 p < 0.05 阈值下通过检验的社区数量。
社区质控结果可用于辅助判断当前社区划分是否合理。如果出现大量规模过小或结构不够紧密的社区，可结合质控结果考虑适当降低 resolution 参数；反之，则可考虑提高该参数，以探索更细粒度的社区结构。参数调整应综合考虑网络拓扑特征、社区规模分布及生物学背景，而不宜仅依据单一指标作出判断。
```python
rdata.community_qc(n_permutations=500)  #质控聚类的社区，观察分辨率会不会太高，导致社区划分过于小，需要调整分辨率参数
```

最后，我们将利用基于 APEX-seq 数据整理的细胞器定位标记基因集，对聚类得到的各个基因共定位社区进行区域富集分析，以评估不同社区的整体细胞器定位倾向性，并探索基因共定位网络的空间功能特征。
```python
gene_loc_file = "/mnt/dfc_data3/project/linyusen/database/17_RNA_colocation/01_marker/data/gene_loc.csv" #细胞器定位标记基因集
rdata.community_enrich_region(gene_loc_file=gene_loc_file) #计算社区富集的细胞区域
```

## 保存结果
完成社区划分后，我们将按照社区逐个导出基因共定位网络，并以 Cytoscape 可读取的文本格式保存。

每个社区输出一个 `community_<社区ID>.txt` 文件，内容为制表符分隔的边表，包含以下字段：
- `source`：边的起始基因
- `target`：边的目标基因
- `p_val`：共定位 p 值
- `effect_area`：效应面积
这些文件可用于后续网络检查或导入 Cytoscape。
```python
output_dir = '/path_to_savedir/'
rdata.save_communities(output_dir=output_dir) #保存社区为cytoscape可导入的文件
```
