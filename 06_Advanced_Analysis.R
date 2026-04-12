# ==========================================
# 第六步：高阶下游分析 (Advanced Analysis)
# 核心目的：解析细胞间的空间邻接关系（微环境）与信号交流（通讯）
# ==========================================

# ------------------------------------------
# 1. 构建空间微环境 (Spatial Niche Analysis)
# 底层逻辑：机器会以每一个细胞为圆心画个圈，看看它周围最常挨着哪些兄弟。
# 经常固定搭配出现的细胞组合，就会被定义为一个“微环境 (Niche)”。
# ------------------------------------------

# 计算细胞邻居关系
# neighbors.k = 30 代表看周围最近的 30 个细胞
# niches.k = 5 代表我们期望初步把切片划分为 5 种微环境（比如：肿瘤核心区、免疫浸润区、间质区等）
xenium.obj <- BuildNicheAssay(object = xenium.obj, 
                              fov = "fov", 
                              group.by = "cell_type", # 基于我们第四步注释好的细胞类型来算
                              niches.k = 5, 
                              neighbors.k = 30)

# 像处理基因表达一样，对“微环境”进行聚类
# 这一步跑完，每个细胞除了拥有“细胞类型”的身份，还多了一个“所属微环境”的标签
xenium.obj <- FindClusters(xenium.obj, assay = "niche", resolution = 0.5)

# 出图验证微环境！
# 你会看到切片被划分成了一块一块的“领地”（比如中间全是红色的一大块肿瘤区，边缘是蓝色的免疫区）
p_niche <- ImageDimPlot(xenium.obj, fov = "fov", group.by = "niche_clusters", axes = TRUE)
print(p_niche)


# ------------------------------------------
# 2. 探究微环境的细胞成分 (Niche Composition)
# 底层逻辑：我们要看看刚才分出来的“0号微环境”里，到底主要是哪些细胞在开会？
# ------------------------------------------
# 画一个柱状图，直观展示每个微环境里不同细胞类型的比例
# 比如你会发现：某一个 Niche 里 90% 都是肿瘤细胞和成纤维细胞，几乎没有 T 细胞，这就是典型的“冷肿瘤”区域。
p_composition <- VlnPlot(xenium.obj, features = "cell_type", group.by = "niche_clusters", pt.size = 0)
# (注：用 Seurat 画成分比例图有很多种进阶代码，这里保留最基础的数据探索抓手)


# ------------------------------------------
# 3. 细胞间通讯的直观证据：受体-配体共定位 (Ligand-Receptor Colocalization)
# 底层逻辑：真正的细胞通讯软件（如 CellChat）代码极其庞大。
# 但在空间组学里，最硬核、最直观的通讯证据，就是证明“发信号的基因(配体)”和“接收信号的基因(受体)”在物理空间上紧紧挨着！
# ------------------------------------------

# 假设我们要研究经典的免疫检查点通路：PD-L1 (基因名 CD274) 和 PD-1 (基因名 PDCD1)
# 我们让 CD274 发红光，PDCD1 发绿光。
# blend = TRUE 是空间通讯最核心的参数！它会把红绿光叠加，如果它们在空间上挨在一起交谈，就会融合成【黄色】！
p_communication <- ImageFeaturePlot(xenium.obj, 
                                    fov = "fov", 
                                    features = c("CD274", "PDCD1"), 
                                    blend = TRUE, 
                                    max.cutoff = "q95",
                                    size = 0.5)

# 这张图一旦跑出来并且有大片的黄色区域，放在文章里就是极其有力的细胞通讯证据
print(p_communication)


# ------------------------------------------
# 4. 批量保存今天的所有成果！
# ------------------------------------------
pdf("06_Niche_and_Communication_Results.pdf", width = 15, height = 8)
print(p_niche)
print(p_communication)
dev.off()

# 【极其重要】：保存你辛辛苦苦处理好的完整 Seurat 对象
# 跑完这六步，这个对象里已经包含了降维、聚类、注释、微环境等所有心血。
# 把它存下来，下次直接 load 就可以从任何一步继续分析，不用从头再跑了！
saveRDS(xenium.obj, file = "Xenium_Fully_Analyzed_Object.rds")