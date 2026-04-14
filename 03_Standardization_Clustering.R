# ==========================================
# 第三步：数据标准化与降维聚类
# 核心目的：消除系统误差，提取核心特征，将细胞分群
# ==========================================

# ------------------------------------------
# 1. 数据标准化 (Normalization)
# ------------------------------------------
# SCTransform 是目前 10x 官方和高分文献最推荐的算法
# 它能一步到位完成数据的标准化、寻找高变基因和数据缩放，极大地保留生物学真实差异
# 但是！！！只适用于数量较少的情况，比如一万细胞左右
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", clip.range = c(-10, 10))

#大体量数据用一下代码
#xenium.obj <- NormalizeData(xenium.obj)
#xenium.obj <- FindVariableFeatures(xenium.obj, selection.method = "vst")
xenium.obj <- ScaleData(xenium.obj)

# ------------------------------------------
# 2. PCA 线性降维 (提取主要特征)
# ------------------------------------------
# 把几千个基因的复杂维度，压缩成最重要的 50 个主成分 (PCs)
# features = rownames(xenium.obj))这个代码把数据对象里**所有的基因（全体名单）**全塞给电脑，让它去算降维
# 可以改用VariableFeatures(object = xenium.obj))，就是只选用关键基因
xenium.obj <- RunPCA(xenium.obj, npcs = 50, VariableFeatures(object = xenium.obj))

# 【关键术间检查】：画出碎石图 (Elbow Plot)
# 运行这行代码后，右下角会出现一张图。图上的点会像一条下垂的手臂。
# 绝大多数情况下，我们需要找到那个“手肘”弯折的位置对应的数字（通常在 15 到 30 之间）。
# 这个数字代表了我们后续分析需要用到的有效维度数量。
ElbowPlot(xenium.obj, ndims = 50)


# ------------------------------------------
# 3. 寻找邻居并聚类 (Clustering)
# ------------------------------------------
# 【注意】：这里的 dims = 1:30 是示例。
# 你需要根据上面那张碎石图（Elbow Plot）的“手肘”位置，来修改这里的 30！
pcs_use <- 1:30 

# 构建细胞邻接图（寻找每个细胞最近的邻居）
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = pcs_use)

# 执行聚类。resolution（分辨率）是这步最核心的参数！
# 绝大多数肿瘤分析中，我们不知道到底分几个群最合适，所以这里直接一次性跑多个分辨率。
# 数值越大，切分得越碎（0.3比较粗糙，0.8分得很细）。
xenium.obj <- FindClusters(xenium.obj, resolution = c(0.3, 0.5, 0.8))

# 设定当前我们要观察的分辨率（比如先看看 0.5 分得怎么样）
Idents(xenium.obj) <- "SCT_snn_res.0.5"


# ------------------------------------------
# 4. UMAP 非线性降维与出图 (Visualization)
# ------------------------------------------
# 把几十维的特征空间，拍扁到人眼能看懂的二维 XY 平面上
xenium.obj <- RunUMAP(xenium.obj, reduction = "pca", dims = pcs_use)

# 术后复查出图！
# 运行后你会看到一张彩色的图，上面标着 0, 1, 2... 的数字。
# 每一个数字代表一个拥有相似基因表达特征的细胞群（可能是肿瘤细胞群、巨噬细胞群等）。
DimPlot(xenium.obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#0.5代表点的大小，可以灵活调整

# 【图 1：经典的 UMAP 调色盘】# pt.size 设小一点，防止黑墙效应,这个代码可能更好
#p1 <- DimPlot(xenium.obj, reduction = "umap", label = TRUE, pt.size = 0.1) +
#NoLegend() + ggtitle("UMAP 细胞分群")# 【图 2：空间原位映射图 (终极杀招！)】# 把 UMAP 上分出来的 0、1、2... 细胞群，按照真实物理坐标“贴”回肿瘤切片上！
#p2 <- ImageDimPlot(xenium.obj, size = 0.5) +
#ggtitle("真实物理空间上的细胞群落分布")# 把两张图并排显示！
