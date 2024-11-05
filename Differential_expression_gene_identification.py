#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Differential_expression_gene_identification.py    
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2023/3/10 上午10:25   skychou      1.0         None
'''
import logging as lg
import os
import click
from rpy2.robjects import r
import pickle


@click.command()
@click.option("--input_file", "-i", required=True, help="输入计数矩阵文件")
@click.option("--output_dir", "-o", default="./DEG_identification_output",
              help="输出质控分析结果目录，默认为'DEG_identification_output'")
@click.option("--software", "-s", type=click.Choice(["edgeR", "DESeq2"]), default="edgeR",
              help="差异分析R包,默认为'edgeR'")
@click.option("--foldchange", "-f", default="1", help="差异分析log2foldchange筛选阈值，默认为1")
@click.option("--padj", "-p", default="0.05", help="差异分析矫正p值padj筛选阈值，默认为0.05")
@click.option("--compare", "-cp", default="", help="需要比较的分组信息（列表的形式，对照组在前，实验组在后）")
@click.option("--config", "-c", type=click.Path(), default="./config.py",
              help="配置文件目录，默认为./config.py")
def main(input_file, output_dir, software, foldchange, padj, compare, config):
    # 获得config文件的内容
    import sys
    global cfg
    if not os.path.exists(config):
        raise Exception(f"配置文件{config}不存在，请检查路径")
    try:
        module_name = os.path.basename(config).split('.')[0]
        cfg = __import__(module_name)
    except:
        sys.path.append(os.path.dirname(config))
        cfg = __import__(module_name)

    global config_prefix
    config_prefix = os.path.splitext(config)[0]
    degi_args = cfg.degi_args
    volcano_color_up = degi_args.get('volcano_color_up')
    volcano_color_down = degi_args.get('volcano_color_down')

    # 设置日志文件
    # 记录器
    global logger
    logger = lg.getLogger("mylog")
    logger.setLevel(lg.DEBUG)

    # 处理器
    consleHandler = lg.StreamHandler()
    consleHandler.setLevel(lg.DEBUG)

    # out 目录检查并处理
    if output_dir == "./":
        output_dir = os.path.dirname(__file__)
    if output_dir != os.path.dirname(__file__):
        if output_dir.endswith("/"):
            output_dir = output_dir[:-1]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    if not os.path.exists(output_dir):
        raise Exception("out dir does not exist!输出路径不存在！")

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    # 没有给定handler日志级别，将用logger的级别
    global log_file
    log_file = f"{output_dir}differential_expression_gene_identification.log"
    fileHandler = lg.FileHandler(filename=log_file, mode="w+", encoding="UTF-8")
    fileHandler.setLevel(lg.INFO)

    # formatter格式
    # formatter = lg.Formatter("%(asctime)s|%(levelname)-8s|%(filename)10s:%(lineno)4s|%(message)s")
    formatter = lg.Formatter("%(asctime)s|%(filename)10s:%(lineno)4s|%(message)s")
    # 给处理器设置格式
    consleHandler.setFormatter(formatter)
    fileHandler.setFormatter(formatter)

    # 给记录器设置处理器
    logger.addHandler(consleHandler)
    logger.addHandler(fileHandler)

    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        logger.error(f"输入文件{input_file}不存在，请检查路径！\n")
        raise Exception(f"输入文件{input_file}不存在，请检查路径！")

    # 检查输出目录是否存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获得分组和样本信息
    import ast
    compares=ast.literal_eval(compare)
    groups = cfg.user_args["groups"]
    groups = {key: groups[key] for key in compares}
    samples = cfg.user_args["samples"]
    d = {}
    for key, value in groups.items():
        for i in value:
            d[i] = key
    templist = []
    templist2 = []
    for key, value in d.items():
        templist.append(f'''"{key}" = "{value}"''')
        templist2.append(f'''"{key}"''')


    gs_vector = f'''gs_vector <- c({",".join(templist)})'''
    choice = f'''c({",".join(templist2)})'''

    if software == "edgeR":
        cmd = f'''
#设置工作目录
work_dir <- "{os.path.dirname(os.path.abspath(__file__))}" #变量
setwd(work_dir)
###edgeR
library(edgeR)
count_file <- "{input_file}" #变量
x <-
  read.table(
    count_file,
    head = TRUE,
    sep = "\\t",
    row.names = 1,
    check.names = FALSE
  )
#修改列名
colnames(x) <-
  unlist(sapply(colnames(x), function(y) sub(".*/(.*)\\\\.genes.*", "\\\\1", y)))
x <- x[,{choice}]
{gs_vector}###新增
tempc <- c()###新增
for (i in names(x)){{tempc <- c(tempc,gs_vector[i])}}###新增
group <- factor(as.vector(tempc)) #这里后面要作为变量，自己识别
group_df <- data.frame(sample = names(x),group = group)
y <- DGEList(counts = x, group = group)
keep <- rowSums(cpm(y) > 1) >= 2 ##至少在两个样本里cpm均大于1
y <- y[keep, , keep.lib.sizes = FALSE]
# Normalization
y <- calcNormFactors(y)
#计算每个基因expected_count值的离散度
#估计变异系数，即估计方差；估计内部差异程度，看组间差异是否比内部差异大，如果大，可选为差异基因
y <- estimateCommonDisp(y, verbose = TRUE)
out_file <- paste(count_file, sep = ".", "normalized")
#out_file <- paste('{output_dir}',out_file,sep = "") ###新增
write.table(y$pseudo.counts, out_file)
print("normalized文件创建成功！")     
pdf("{output_dir}plotMDS2.pdf")#变量
par(
  mfrow = c(2, 2),
  cex = 0.5,
  cex.axis = 1,
  cex.lab = 1,
  cex.main = 1.2
)
plotMDS(y, main = "method:logFC", col = "red")
plotMDS(y,
        method = "bcv",
        main = "method:bcv",
        col = "red")
dev.off()
print("MDS文件创建成功！")

## 聚类分析
data = t(y$pseudo.counts)
data2 = scale(t(y$pseudo.counts), center = T, scale = T) #数据标准化,除以方差
#层次聚类
d <- dist(data, method = "euclidean")
d2 <- dist(data2, method = "euclidean")
hc <- hclust(d, "single")
hc2 <- hclust(d2, "single")
pdf("{output_dir}clustering_tree_plot.pdf")#变量
par(
  mfrow = c(2, 2),
  cex = 1,
  cex.axis = 1,
  cex.lab = 1,
  cex.main = 1.2
)
plot(hc, xlab = "samples")
plot(hc2, main = "Cluster Dendrogram(scale)", xlab = "samples")
dev.off()
print("样本聚类分析完成！")


y <- estimateTagwiseDisp(y)

#获得样本间的相似性数据
c <- cor(y$counts) #计算样本间的相关系数矩阵

##差异分析
#进行精确检验=》以此来筛选差异表达基因DEGs
et <- exactTest(y)

#just all genes
ordered_tags <- topTags(et, n = length(et$table[, 1]))
#extract DEGs
allDEGs = ordered_tags$table
#remove lines of calculation errors of FDR value
allDEGs = allDEGs[is.na(allDEGs$FDR) == FALSE, ]
length(allDEGs[, 1]) #[1] 14831
out_file2 = "{output_dir}edgeR_allGenes.out" #变量
write.table(allDEGs, out_file2)

#设定阈值
foldchange = {foldchange} #logFC =1, which is equivalent to 2 times
padj = {padj} #FDR = 0.05

#significant DEGs
sigDEGs = allDEGs[allDEGs$FDR < padj, ]
#significant DEGs and logFC >1
sigDEGs2 = sigDEGs[(sigDEGs$logFC > foldchange |
                      sigDEGs$logFC < (-foldchange)), ]
#2 upregulated times
upDEGs = sigDEGs[sigDEGs$logFC > foldchange, ]
#2 downregulated times
downDEGs = sigDEGs[sigDEGs$logFC < (-foldchange), ]
#output the results of DEGs
out_file3 = "{output_dir}edgeR_sigDEGs_FDR{padj}_logFC{foldchange}.out"
write.table(sigDEGs2, out_file3)
out_file4 = "{output_dir}edgeR_upDEGs_FDR{padj}_logFC{foldchange}.out"
write.table(upDEGs, out_file4)
out_file5 = "{output_dir}edgeR_downDEGs_FDR{padj}_logFC{foldchange}.out"
write.table(downDEGs, out_file5)
out_file6 = "{output_dir}edgeR_sigDEGs_FDR{padj}.out"
write.table(sigDEGs, out_file6)
print("差异分析.out文件生成完成！")
write.table(data.frame(query_id=row.names(sigDEGs),logFC=sigDEGs$logFC),"{output_dir}sigDEG.tsv",row.names = F,quote = FALSE,sep = "\\t")

#绘制火山图
#draw the volcano map of significant DEGs
pdf(file = "{output_dir}volcano.pdf")#变量
xMax = max(allDEGs$logFC) + 1
yMax = 20
plot(
  allDEGs$logFC,
  -log10(allDEGs$FDR),
  xlab = "logFC",
  ylab = "-log10(FDR)",
  main = "Volcano",
  xlim = c(-xMax, xMax),
  ylim = c(0, yMax),
  pch = 20,
  cex = 1
)
points(
  upDEGs$logFC,
  -log10(upDEGs$FDR),
  pch = 20,
  col = "{volcano_color_up}",
  cex = 1
)
points(
  downDEGs$logFC,
  -log10(downDEGs$FDR),
  pch = 20,
  col = "{volcano_color_down}",
  cex = 1
)
abline(v = 0, lty = 2, lwd = 2)
dev.off()
print("火山图绘制完成！")

#绘制散点图
# Drawing of scatter map
#extract normalized expression levels of all genes
norm_ge = y$pseudo.counts
#draw scatter map after log conversion of normalized value
pdf(file = "{output_dir}scatter.pdf")#变量
#normalized expression levels of all genes in Con group
samples1 = group_df$sample[group_df$group == levels(group)[1]]
x_ge = 0
for (i in samples1) {{
  x_ge = x_ge + norm_ge[,i]
}}
x_ge = x_ge / 2 + 1
# x_ge = (norm_ge[, 1] + norm_ge[, 2]) / 2 + 1  #:防止0值对log的影响

#normalized expression levels of all genes in Sam group
samples2 = group_df$sample[group_df$group == levels(group)[2]]
y_ge = 0
for (i in samples2) {{
  y_ge = y_ge + norm_ge[,i]
}}

y_ge = y_ge / 2 + 1
# y_ge = (norm_ge[, 3] + norm_ge[, 4]) / 2 + 1 #防止0值对log的影响

#calculate the ranges of x and y axises
xMin = min(log(x_ge), log(y_ge)) - 1
xMax = max(log(x_ge), log(y_ge)) + 1
yMin = min(log(y_ge), log(y_ge)) - 1
yMax = max(log(y_ge), log(y_ge)) + 1
#plot  the points of expression levels of all genes
plot(
  log(x_ge),
  log(y_ge),
  xlab = levels(group)[1],
  ylab = levels(group)[2],
  main = "Scatter",
  xlim = c(xMin, xMax),
  ylim = c(yMin, yMax),
  pch = 20,
  cex = 1,
  cex.axis = 1.5,
  cex.lab = 1.5
)
#all significant DEGs and logFC >1
#normalized expression levels of upregulated genes in Con group
up_ge_x = x_ge[rownames(upDEGs)]
#normalized expression levels of upregulated genes in Sam group
up_ge_y = y_ge[rownames(upDEGs)]
#normalized expression levels of downregulated genes in Con group
down_ge_x = x_ge[rownames(downDEGs)]
#normalized expression levels of downregulated genes in Sam group
down_ge_y = y_ge[rownames(downDEGs)]
#plot  the points of significant DEGs
points(
  log(up_ge_x),
  log(up_ge_y),
  pch = 20,
  col = "{volcano_color_up}",
  cex = 1
)
points(
  log(down_ge_x) ,
  log(down_ge_y),
  pch = 20,
  col = "{volcano_color_down}",
  cex = 1
)
abline(a = 0,
       b = 1,
       lty = 2,
       lwd = 1) #append a line of x_ge = y_ge
dev.off()
print("散点图绘制完成！")

#聚类分析热图
#extract normalized data of significant DEGs
hmData = norm_ge[rownames(sigDEGs), ]
length(hmData[, 1]) #[1] 211=> logFC>1(in this case)
hmMat = as.matrix(hmData)
library(pheatmap)
pdf(file = "{output_dir}sigDEGs_pheatmap.pdf")
pheatmap(
  hmMat,
  scale = "row",
  clustering_distance_row = "correlation",
  fontsize = 10,
  fontsize_row = 1
)
dev.off()

print("聚类分析热图绘制完成！")
'''
    if software == "DESeq2":
        cmd = f'''
work_dir <- "{os.path.dirname(os.path.abspath(__file__))}" #变量
setwd(work_dir)
##DESeq2
library(DESeq2)
library(limma)
count_file <- "{input_file}" #变量
x <-
  read.table(
    count_file,
    head = TRUE,
    sep = "\\t",
    row.names = 1,
    check.names = FALSE
  )

#修改列名
colnames(x) <-
  unlist(sapply(colnames(x), function(y) sub(".*/(.*)\\\\.genes.*", "\\\\1", y)))
x <- x[,{choice}]
group <- factor(c("group1", "group1", "group2", "group2")) #这里后面要作为变量，自己识别
{gs_vector}###新增
tempc <- c()###新增
for (i in names(x)){{tempc <- c(tempc,gs_vector[i])}}###新增
group <- factor(as.vector(tempc)) #这里后面要作为变量，自己识别
group_df <- data.frame(sample = names(x),group = group)
coldata <- data.frame(condition = factor(rep(c("Con", "Sam"), each = 2), levels = c('Con', 'Sam')))
# 构建DESeqDataSet对象
x <- round(x)
dds <-
  DESeqDataSetFromMatrix(
    countData = x,
    colData = coldata,
    design = ~ condition
  )

#进行归一化和差异分析
dds <- DESeq(dds)

# 过滤低表达量的基因
keep <- rowSums(counts(dds) >= 10) >= 2 #至少在2个样本中count大于等于10
dds <- dds[keep,]

# 绘制波利斯图（横轴是对数倍数变化，纵轴是方差/均值比例）
# plotDispEsts(dds)

# 获取标准化后的表达值
norm_counts <- counts(dds, normalized = TRUE)
out_file <- paste(count_file, sep = ".", "normalized") #变量
#out_file <- paste('{output_dir}',out_file,sep = "") ###新增
write.table(norm_counts, out_file)
print("normalized文件创建成功！")
#MDS图
pdf("{output_dir}plotMDS2.pdf")#变量
par(
  mfrow = c(2, 2),
  cex = 0.5,
  cex.axis = 1,
  cex.lab = 1,
  cex.main = 1.2
)
plotMDS(dds,main = "method:logFC", col = "red")
plotMDS(dds,
        method = "bcv",
        main = "method:bcv",
        col = "red")
dev.off()
print("MDS文件创建成功！")

## 聚类分析
data = t(norm_counts)
data2 = scale(t(norm_counts),center = TRUE,scale = TRUE) #数据标准化,除以方差
d <- dist(data, method = "euclidean")
d2 <- dist(data2, method = "euclidean")
hc <- hclust(d, "single")
hc2 <- hclust(d2, "single")
pdf("{output_dir}clustering_tree_plot.pdf")#变量
par(
  mfrow = c(2, 2),
  cex = 1,
  cex.axis = 1,
  cex.lab = 1,
  cex.main = 1.2
)
plot(hc, xlab = "samples")
plot(hc2, main = "Cluster Dendrogram(scale)", xlab = "samples")
dev.off()

#获得样本间的相似性数据
c <- cor(norm_counts) #计算样本间的相关系数矩阵

# 获取差异表达基因结果
allDEGs <- results(dds)
allDEGs <-  as.data.frame(allDEGs)
allDEGs = allDEGs[is.na(allDEGs$padj) == FALSE,]

foldchange = 1 #logFC =1, which is equivalent to 2 times
padj = 0.05 #FDR = 0.05

# 筛选显著差异表达基因
sigDEGs <- allDEGs[which(allDEGs$padj < padj), ]
sigDEGs

#significant DEGs and logFC >1
sigDEGs2 = sigDEGs[(sigDEGs$log2FoldChange > foldchange |
                      sigDEGs$log2FoldChange < (-foldchange)),]

# 筛选上调或下调基因
upDEGs <- sigDEGs[which(sigDEGs$log2FoldChange > foldchange), ]
downDEGs <- sigDEGs[which(sigDEGs$log2FoldChange < -foldchange), ]

# 输出差异表达基因结果
out_file2 = "{output_dir}DESeq2_allGenes.out" #变量
write.table(allDEGs, out_file2)
out_file3 = "{output_dir}DESeq2_sigDEGs_FDR0.05_logFC1.out" #变量
write.table(sigDEGs2, out_file3)
out_file4 = "{output_dir}DESeq2_upDEGs_FDR0.05_logFC1.out" #变量
write.table(upDEGs, out_file4)
out_file5 = "{output_dir}DESeq2_downDEGs_FDR0.05_logFC1.out" #变量
write.table(downDEGs, out_file5)
out_file6 = "{output_dir}DESeq2_sigDEGs_FDR0.05.out"
write.table(sigDEGs, out_file6)

#保存差异基因列表
write.table(data.frame(query_id=row.names(sigDEGs2),logFC=sigDEGs2$log2FoldChange),"{output_dir}sigDEG.tsv",row.names = F,quote = FALSE,,sep = "\\t")

#绘制火山图
pdf(file = "{input_file}volcano.pdf") #变量
xMax = max(allDEGs$log2FoldChange) + 1
yMax = 20
plot(
  allDEGs$log2FoldChange,-log10(allDEGs$padj),
  xlab = "logFC",
  ylab = "-log10(FDR)",
  main = "Volcano",
  xlim = c(-xMax, xMax),
  ylim = c(0, yMax),
  pch = 20,
  cex = 1
)
points(
  upDEGs$log2FoldChange,-log10(upDEGs$padj),
  pch = 20,
  col = "{volcano_color_up}",
  cex = 1
)
points(
  downDEGs$log2FoldChange,-log10(downDEGs$padj),
  pch = 20,
  col = "{volcano_color_down}",
  cex = 1
)
abline(v = 0, lty = 2, lwd = 2)
dev.off()

#绘制散点图
# Drawing of scatter map
#extract normalized expression levels of all genes
norm_ge = norm_counts
#draw scatter map after log conversion of normalized value
pdf(file = "{input_file}scatter.pdf")#变量
#normalized expression levels of all genes in Con group
samples1 = group_df$sample[group_df$group == levels(group)[1]]
x_ge = 0
for (i in samples1) {{
  x_ge = x_ge + norm_ge[,i]
}}
x_ge = x_ge / 2 + 1
# x_ge = (norm_ge[, 1] + norm_ge[, 2]) / 2 + 1  #:防止0值对log的影响

#normalized expression levels of all genes in Sam group
samples2 = group_df$sample[group_df$group == levels(group)[2]]
y_ge = 0
for (i in samples2) {{
  y_ge = y_ge + norm_ge[,i]
}}
y_ge = y_ge / 2 + 1
# y_ge = (norm_ge[, 3] + norm_ge[, 4]) / 2 + 1 #防止0值对log的影响
#calculate the ranges of x and y axises
xMin = min(log(x_ge), log(y_ge)) - 1
xMax = max(log(x_ge), log(y_ge)) + 1
yMin = min(log(y_ge), log(y_ge)) - 1
yMax = max(log(y_ge), log(y_ge)) + 1
#plot  the points of expression levels of all genes
plot(
  log(x_ge),
  log(y_ge),
  xlab = "Con",
  ylab = "Sam",
  main = "Scatter",
  xlim = c(xMin, xMax),
  ylim = c(yMin, yMax),
  pch = 20,
  cex = 1,
  cex.axis = 1.5,
  cex.lab = 1.5
)
#all significant DEGs and logFC >1
#normalized expression levels of upregulated genes in Con group
up_ge_x = x_ge[rownames(upDEGs)]
#normalized expression levels of upregulated genes in Sam group
up_ge_y = y_ge[rownames(upDEGs)]
#normalized expression levels of downregulated genes in Con group
down_ge_x = x_ge[rownames(downDEGs)]
#normalized expression levels of downregulated genes in Sam group
down_ge_y = y_ge[rownames(downDEGs)]
#plot  the points of significant DEGs
points(
  log(up_ge_x),
  log(up_ge_y),
  pch = 20,
  col = "{volcano_color_up}",
  cex = 1
)
points(
  log(down_ge_x) ,
  log(down_ge_y),
  pch = 20,
  col = "{volcano_color_down}",
  cex = 1
)
abline(a = 0,
       b = 1,
       lty = 2,
       lwd = 1) #append a line of x_ge = y_ge
dev.off()

#聚类分析热图
#extract normalized data of significant DEGs
hmData = norm_ge[rownames(sigDEGs), ]
length(hmData[, 1]) #[1] 211=> logFC>1(in this case)
hmMat = as.matrix(hmData)
library(pheatmap)
pdf(file = "{input_file}sigDEGs_pheatmap.pdf")
pheatmap(
  hmMat,
  scale = "row",
  clustering_distance_row = "correlation",
  fontsize = 10,
  fontsize_row = 1
)
dev.off()
        '''
    try:
        r(cmd)
        logger.info(f"差异分析完成！结果目录为：{output_dir}\n")
    except Exception as e:
        logger.error(f"差异分析失败！报错信息如下：\n{e}")

    # 读取原始 pkl 文件中的列表数据
    config_pkl = f'{config_prefix}.pkl'
    with open(config_pkl, "rb") as f:
        data = pickle.load(f)

    # 在原始列表中添加新的字符串元素
    # 获取out下的所有文件和文件夹
    all_items = os.listdir(output_dir)
    # 过滤掉文件夹，只保留所有的文件
    all_items = [os.path.join(output_dir, i) for i in all_items]
    data["deg_identification_files"] = all_items

    # 将原始列表和新的字符串一起写入新的 pkl 文件中
    with open(config_pkl, "wb") as f:
        pickle.dump(data, f)


if __name__ == '__main__':
    main()
