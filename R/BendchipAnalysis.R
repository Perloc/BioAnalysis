# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' 读取芯片数据并进行标准化和背景校正
#'
#' @param filepath 芯片数据路径
#'
#' @return 返回表达矩阵
#' @export
#'
#' @examples
read.bendchip <- function(filepath) {
  expr_set <- read.ilmn(datafile)
  # 标准化 + 背景校正
  expr_set <- neqc(expr_set)
  dat <- expr_set$E
  dat <- as.data.frame(dat)

  return(dat)
}



#' 去除表达矩阵中重复的基因
#'
#' 保留表达值平均数最大的基因
#'
#' @param dat 最后两列为Symbol和ENTREZID的表达矩阵
#'
#' @return 去重并且以Symbol为行名的表达矩阵
#' @export
#'
#' @examples
removeDuplication <- function(dat) {
  ### 处理一个gene名对应多个探针情况，平均数排序取最大
  dat$mean <- apply(dat[, -c(ncol(dat), ncol(dat)-1)], 1, mean)
  dat <- dat[order(dat$Symbol, dat$mean, decreasing = T),]
  # 按symbol取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果
  dat <- dat[!duplicated(dat$Symbol),]
  dat <- dat[order(dat$ENTREZID, dat$mean, decreasing = T),]
  dat <- dat[!duplicated(dat$ENTREZID),]
  # 把symbol这一列中的每一行给dat作为dat的行名
  rownames(dat) <- dat$Symbol
  dat <- dat[, -c(ncol(dat), ncol(dat)-2)]
  return(dat)
}

#' 画出聚类图
#'
#' @param expr_mat 表达矩阵
#' @param clust_method 聚类方法，可选 "complete", "single", "average", "ward.D"等
#'
#' @export
#'
#' @examples
clust <- function(expr_mat, clust_method = "complete") {
  dat <- t(expr_mat)
  dat <- as.data.frame(dat)
  datExpr_tree <- hclust(dist(dat), method = clust_method)
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub = "", xlab = "",
       cex.axis = 0.9, cex.main = 1.2, cex.lab = 1, cex = 0.7)
}

computeSoftThreshold <- function(dat) {
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 5)
  # 画图
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 <- 0.9;
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
       type="n", main = paste("Scale independence"));
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       labels=powers, cex=cex1, col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity",
       type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1,col="red")
  power <- sft$powerEstimate

  # 样本数小于20，设置为9；20～30，设置为8；30以上，设置为7
  if (is.na(power)) {
    if (nrow(dat) < 20) { power <- 9 }
    else if (nrow(dat) < 30) { power <- 8 }
    else { power <- 7 }
  }

  return(power)
}

#' WGCNA分析
#'
#' @param expr_mat 表达矩阵
#' @param traits_mat 性状矩阵
#'
#' @return 分析结果
#' @export
#'
#' @examples
WGCNA <- function(expr_mat, traits_mat) {
  # 开启多线程
  allowWGCNAThreads()

  dat <- t(expr_mat)
  power <- computeSoftThreshold(dat)

  cor <- WGCNA::cor
  net = blockwiseModules(
    dat,
    power = power,
    maxBlockSize = 5000,
    TOMType = "unsigned", minModuleSize = 200,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    verbose = 3
  )

  # 可视化
  moduleColors = labels2colors(net$colors)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)

  MEs0 <- moduleEigengenes(dat, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")

  geneTraitCor <- as.data.frame(cor(dat, traits_mat, use = "p"))

  return(list(moduleColors = moduleColors,
              geneTraitCor = geneTraitCor,
              power = power))
}

#' Title
#'
#' 根据cut高度合并module
#'
#' @param expr_mat 表达矩阵
#' @param WGCNA_result WGCNA分析结果
#' @param cutHeightThres cut阈值
#'
#' @return
#' @export
#'
#' @examples
WGCNAMergeModule <- function(expr_mat, WGCNA_result, cutHeightThres = 0.4) {
  dat <- t(expr_mat)
  moduleColors <- WGCNA_result$moduleColors

  merge <- mergeCloseModules(dat, moduleColors, cutHeight = cutHeightThres, verbose = 3)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs

  geneModuleMembership <- as.data.frame(cor(dat, mergedMEs, use = "p"))

  WGCNA_result$geneModuleMembership <- geneModuleMembership
  WGCNA_result$moduleColors <- mergedColors
  WGCNA_result$MEs <- mergedMEs
  return(WGCNA_result)
}

#' Title
#'
#' 画出module与性状的相关性热图
#'
#' @param WGCNA_result WGCNA分析结果
#' @param traits_mat 性状矩阵
#'
#' @return
#' @export
#'
#' @examples
WGCNAHeatmap <- function(WGCNA_result, traits_mat) {
  MEs <- WGCNA_result$MEs
  moduleTraitCor <- WGCNA_result$geneTraitCor

  moduleTraitCor <- cor(MEs, traits_mat, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits_mat))
  sizeGrWindow(10, 6)
  # Will display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  # 模块与性状的相关性热图
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(traits_mat),
                 yLabels = row.names(moduleTraitCor),
                 ySymbols = row.names(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}

#' 画WGCNA分析的散点图
#'
#' @param result WGCNA分析结果
#' @param traits_mat 性状矩阵
#' @param module 模块颜色
#' @param pheno 性状名
#'
#' @return
#' @export
#'
#' @examples
WGCNAScatterPlot <- function(result, traits_mat, module, pheno) {
  modNames <- substring(colnames(result$MEs), 3)
  module_column <- match(module, modNames)
  pheno_column <- match(pheno, colnames(traits_mat))
  moduleGenes <- result$moduleColors == module

  # 画图
  sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  verboseScatterplot(abs(result$geneModuleMembership[moduleGenes, module_column]),
                     abs(result$geneTraitCor[moduleGenes, pheno_column]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", pheno),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

#' 获取指定模块内的基因
#'
#' @param expr_mat 表达矩阵
#' @param result WGCNA分析结果
#' @param module 模块颜色
#'
#' @return
#' @export
#'
#' @examples
getModuleGenes <- function(expr_mat, result, module) {
  probes <- rownames(expr_mat)
  inModule <- (result$moduleColors %in% module)
  modProbes <- probes[inModule]
  return(expr_mat[rownames(expr_mat) %in% modProbes, ])
}

### 差异分析 ###
#' limma差异分析
#'
#' @param expr_mat 表达矩阵
#' @param traits_mat 性状矩阵
#' @param contrasts 对比项
#' @param item 获取哪个对比项的差异分析结果
#' @param pcutoff pvalue阈值
#' @param lfc log foldchange阈值
#' @param title 火山图标题
#'
#' @return
#' @export
#'
#' @examples
limmaAnalysis <- function(expr_mat, traits_mat, contrasts, item = 1, pcutoff = 0.05, lfc = 0,
                          title = "Volcano Plot") {
  fit <- lmFit(expr_mat, traits_mat)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2, trend=TRUE)

  summary(decideTests(fit2, method = "global", adjust.method = "none",
                      p.value = pcutoff, lfc = lfc))

  result <- topTable(fit2, coef = item, number = Inf, sort.by = "p")
  result$change <- ifelse(result$P.Value > pcutoff,'stable',
                          ifelse( result$logFC > lfc,'up',
                                  ifelse( result$logFC < -lfc,'down','stable')))

  result <- result %>% mutate(genelabels = "")
  up <- subset(result, change=="up")
  up$genelabels[1:5] <- as.character(rownames(up)[1:5])
  down <- subset(result, change=="down")
  down$genelabels[1:5] <- as.character(rownames(down)[1:5])

  p <- ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = change)) +
    ggtitle(title) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
    geom_vline(xintercept = c(-1, 1),lty = 4,col = "black",lwd = 0.8)+ #logFC分界线
    geom_hline(yintercept=-log10(0.05),lty = 4,col = "black",lwd = 0.8)+ #adj.p.val分界线
    theme_bw()
  p <- p + geom_text_repel(up, mapping = aes(x = logFC, y = -log10(P.Value), label = genelabels))
  p <- p + geom_text_repel(down, mapping = aes(x = logFC, y = -log10(P.Value), label = genelabels))
  ggsave("Volcano_Plot.pdf", p)

  return(list(result = result, up = up, down = down))
}

#' 根据基因名匹配EntrezID
#'
#' @param DA_result limmaAnalysis返回的结果
#' @param anno 标注信息，包含Symbol和ENTREZID列
#'
#' @return
#' @export
#'
#' @examples
matchEntrezID <- function(DA_result, anno) {
  result <- DA_result$result
  up <- DA_result$up
  down <- DA_result$down

  result$Symbol <- row.names(result)
  result <- inner_join(result, anno, by = "Symbol", multiple = "first")
  up$Symbol <- row.names(up)
  up <- inner_join(up, anno, by = "Symbol", multiple = "first")
  down$Symbol <- row.names(down)
  down <- inner_join(down, anno, by = "Symbol", multiple = "first")

  return(list(result = result, up = up, down = down))
}

#' GSEA分析
#'
#' @param DA_result limma差异分析结果
#'
#' @return
#' @export
#'
#' @examples
GSEAnalysis <- function(DA_result) {
  foldchanges <- DA_result$result$logFC
  names(foldchanges) <- DA_result$result$ENTREZID
  foldchanges <- sort(foldchanges, decreasing = TRUE)

  set.seed(123456)
  gseaKEGG <- gseKEGG(geneList = foldchanges,
                      organism = "hsa",
                      minGSSize = 20,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)

  dotplot(gseaKEGG, split = ".sign") + facet_grid(~.sign)

  # 提取GSEA的结果
  gseaKEGG_results <- gseaKEGG@result
  # 将结果写入csv文件
  write.csv(gseaKEGG_results, "GSEAresults.csv", quote=F)
}
