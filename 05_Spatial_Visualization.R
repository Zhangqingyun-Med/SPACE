# ==========================================
# 第五步：空间可视化 (Spatial Visualization)
# 核心目的：将细胞分类和基因表达映射回二维物理组织切片上
# ==========================================

# 确保使用的身份标签是我们第四步刚注释好的真实细胞名字
Idents(xenium.obj) <- "cell_type"

# ------------------------------------------
# 1. 战局全景：绘制包含所有细胞类型的空间分布图
# ------------------------------------------
# ImageDimPlot 是空间组学最核心的画图函数
# fov = "fov" 代表视场（Field of View），Xenium 数据默认叫这个名字
# axes = TRUE 会在图片边缘画出物理坐标轴，方便评估切片大小
# cols = "polychrome" 是一种适合展示多种细胞的鲜艳配色方案
p1 <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE, cols = "polychrome", size = 0.1)

# 在 RStudio 右下角查看全景图
print(p1)


# ------------------------------------------
# 2. 精确打击：只高亮观察特定的细胞亚群 (Highlight)
# ------------------------------------------
# 在极其密集的肿瘤切片里，细胞全都挤在一起，全景图往往看着很乱。
# 这时候我们需要把无关细胞变成灰色背景，只点亮我们关心的几个群体。
# 比如：我们想看看 T 细胞（杀手）和肿瘤细胞（靶子）的空间位置关系
#在这里加上一个命名代码，一定核对名称是否一样
tumor_cells <- WhichCells(xenium.obj, idents = c("Tumor cells (EPCAM+)", "Tumor cells (LAPTM4B+)", "Proliferating Tumor"))
t_cells <- WhichCells(xenium.obj, idents = c("T cells", "T cells (CXCR4+)"))

highlight_list <- list(
  "Cancer" = tumor_cells, 
  "T_Immune" = t_cells
)
p2 <- ImageDimPlot(xenium.obj, 
                   fov = "fov", 
                   cells.highlight = highlight_list, 
                   cols.highlight = c("red", "blue"), # 对应上面列表里的两个组
                   size = 0.1) +                      # 22万细胞，点一定要小
      ggtitle("肿瘤实质与 T 细胞空间交锋图")
print(p2)


# ------------------------------------------
# 3. 靶点显影：观察特定基因在切片上的表达强度分布
# ------------------------------------------
# ImageFeaturePlot 用来画基因（而不是细胞群）。
# 颜色越亮（通常是红/黄），代表这个基因在这个物理位置表达量越高。
# 比如：我们看看缺氧标志物(HIF1A)或者某个免疫检查点(CD274/PD-L1)在哪里亮起
features_to_plot <- c("EPCAM", "CD8A") # 替换成你感兴趣的基因名字

p3 <- ImageFeaturePlot(xenium.obj, 
                       fov = "fov", 
                       features = features_to_plot, 
                       max.cutoff = "q95", # 屏蔽极少数高得离谱的噪点，让图片对比度更好
                       size = 0.1)
print(p3)


# ------------------------------------------
# 4. 【高阶实操】：保存适合发文章的高清大图
# ------------------------------------------
# Xenium 的切片非常大，动辄几十万个细胞，RStudio 自带的画图窗口往往会被挤得模糊不清。
# 课题组汇报时，必须用代码将其导出为大尺寸的 PDF 或高分辨率 PNG。

# 导出一张宽10英寸、高8英寸的高清 PDF（科研界最喜欢 PDF，因为无限放大不失真）
pdf("Spatial_AllCells_Map.pdf", width = 10, height = 8)
print(p1)
dev.off() # 这行代码极其关键，意思是“关闭画板，保存文件”

# 如果你想给师兄微信发图，也可以导出一张高清 PNG
png("Spatial_Marker_Expression.png", width = 2000, height = 1500, res = 300)
print(p3)
dev.off()
