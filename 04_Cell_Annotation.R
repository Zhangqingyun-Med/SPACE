# ==========================================
# 第四步：细胞类型注释 (Cell Type Annotation)
# 核心目的：给冷冰冰的数字分群赋予真实的生物学身份
# ==========================================

# ------------------------------------------
# 1. 寻找每个群的“专属身份证” (Marker 基因)
# ------------------------------------------
# FindAllMarkers 会计算出每个群里表达量显著高于其他群的基因
# 【注意】：这一步计算量比较大，如果细胞多可能需要跑几分钟，RStudio 出现红灯停顿是正常的，耐心等待。
# only.pos = TRUE: 只找高表达的基因（我们通常不关心哪个基因表达低）
# min.pct = 0.25: 这个基因至少要在该群 25% 的细胞中表达
# logfc.threshold = 0.25: 表达差异倍数阈值
# 加上 max.cells.per.ident = 1000，电脑从每个群里最多随机抽取 1000 个精锐去对比，极大节省时间，且完全不影响找核心 Marker！
all.markers <- FindAllMarkers(xenium.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25，max.cells.per.ident = 1000)

# 【极其推荐的高阶操作】：把找出来的 Marker 基因表保存到本地电脑！
# 这样你可以用 Excel 打开它，用平时学的医学知识，慢慢看每个群最高表达的是什么基因。
write.csv(all.markers, file = "Cluster_Markers_Result.csv")


# ------------------------------------------
# 2. 查阅文献与常识，人工鉴定细胞身份 (Manual Annotation)
# ------------------------------------------
# 假设你看完刚刚导出的 Excel 表，结合经典靶点发现：
# 0群高表达 EPCAM, KRT18 -> 明确是肿瘤细胞/上皮细胞 (Tumor/Epithelial)
# 1群高表达 CD3D, CD8A  -> 明确是 T 细胞 (T cells)
# 2群高表达 CD68, CD163 -> 明确是巨噬细胞 (Macrophages)
# 3群高表达 COL1A1, FAP -> 明确是成纤维细胞 (Fibroblasts)
# 4群高表达 PECAM1, VWF -> 明确是内皮细胞 (Endothelial)

# 构建一个“字典”，把数字映射成真实名字 
# 【注意】：这里的顺序必须和你的群号（0, 1, 2...）完全对应！
# 下面的 5 个群仅仅是通用演示，真实跑出来可能有 10 个甚至 15 个群，你需要相应地加上去。
cell_type_dict <- c(
  "0" = "Tumor cells",
  "1" = "T cells",
  "2" = "Macrophages",
  "3" = "Fibroblasts",
  "4" = "Endothelial"
)


# ------------------------------------------
# 3. 执行“偷梁换柱”，正式赋予身份
# ------------------------------------------
xenium.obj <- RenameIdents(xenium.obj, cell_type_dict)

# 把改好的名字存进对象的元数据表（也就是你文档里提到的“细胞属性表”）里，方便以后随时调用
xenium.obj$cell_type <- Idents(xenium.obj)


# ------------------------------------------
# 4. 术后出图：带有真实细胞名字的二维分布与验证
# ------------------------------------------
# 图 1：画出带有名字的 UMAP 图，现在图上的标签变成了 Tumor, T cells 等
DimPlot(xenium.obj, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4) + NoLegend()+ ggtitle("UMAP 细胞鉴定图")
# 图 2：带名字的空间切片 (Xenium 的灵魂所在)
p2 <- ImageDimPlot(xenium.obj, size = 0.5) + 
      NoLegend() + ggtitle("切片原位细胞分布")
# 并排震撼出图
p1 + p2

# 【高阶验证图】：气泡图 (DotPlot)
# 拿你最关心的几个经典基因去验证一下，看看是不是真的只在对应的细胞群里高表达
# 圆圈越大代表表达该基因的细胞比例越高，颜色越红代表表达量越高
features_to_check <- c("EPCAM", "CD3D", "CD68", "COL1A1", "PECAM1")
DotPlot(xenium.obj, features = features_to_check) + RotatedAxis()
