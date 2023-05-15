library(clusterProfiler)
library(ggplot2)
library(ggthemes)
library(org.Hs.eg.db)


go_kegg_plot <- function(up.data, down.data, method = "GO", n = 5) {
  if (n < nrow(up.data)) {
    up.data <- up.data[1:n, ]
  }
  if (n < nrow(down.data)) {
    down.data <- down.data[1:n, ]
  }

  up.data$change <- "up"
  down.data$change <- "down"
  updown <- rbind(up.data,down.data)
  updown$pl <- ifelse(updown$change == "up",-log10(updown$p.adjust),log10(updown$p.adjust))

  updown <- updown[order(updown$pl,decreasing = F),]
  updown$Description = factor(updown$Description,levels = unique(updown$Description),ordered = TRUE)

  gk_plot <- ggplot(updown,aes(reorder(Description, pl), y=pl)) +
    geom_bar(aes(fill=factor((pl>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
    coord_flip() +
    scale_fill_manual(values=c("#0072B2", "#B20072"), guide=FALSE) +
    labs(x="", y="" ) +
    theme_pander()  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.line.x = element_line(linewidth = 0.3, colour = "black"),#x轴连线
          axis.ticks.length.x = unit(-0.20, "cm"),#修改x轴刻度的高度，负号表示向上
          axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),##线与数字不要重叠
          axis.ticks.x = element_line(colour = "black",linewidth = 0.3) ,#修改x轴刻度的线
          axis.ticks.y = element_blank(),
          axis.text.y  = element_text(hjust=0),
          panel.background = element_rect(fill=NULL, colour = 'white')
    )
  ggsave(paste(method, "_histplot.pdf", sep = ""), gk_plot)
}


#' 对up和down分别进行GO/KEGG分析
#'
#' @param up 带EntrezID列的上调基因
#' @param down 带EntrezID列的下调基因
#' @param method "GO" / "KEGG"
#' @param n 图中显示项的数量
#'
#' @return
#' @export
#'
#' @examples
enrich_split <- function(up, down, method = "GO", n = 5) {
  if (method == "GO") {
    cluster_up <- enrichGO(gene = up$ENTREZID,
                           # universe = all_genes,
                           keyType = "ENTREZID", # ENSEMBL, SYMBOL
                           OrgDb = org.Hs.eg.db,
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
    cluster_down <- enrichGO(gene = down$ENTREZID,
                             # universe = all_genes,
                             keyType = "ENTREZID", # ENSEMBL, SYMBOL
                             OrgDb = org.Hs.eg.db,
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             readable = TRUE)
  } else {
    cluster_up <- enrichKEGG(gene = up$ENTREZID,
                             use_internal_data = F,
                             keyType = 'kegg',  # KEGG 富集
                             organism = 'hsa',  # 物种名称
                             pAdjustMethod = 'none',  # 指定p值校正方法
                             pvalueCutoff = 0.05,  #指定p值阈值（可指定 1 以输出全部）
                             qvalueCutoff = 0.2)  #指定q值阈值（可指定 1 以输出全部）
    cluster_down <- enrichKEGG(gene = down$ENTREZID,
                               use_internal_data = F,
                               keyType = 'kegg',  # KEGG 富集
                               organism = 'hsa',  # 物种名称
                               pAdjustMethod = 'none',  # 指定p值校正方法
                               pvalueCutoff = 0.05,  #指定p值阈值（可指定 1 以输出全部）
                               qvalueCutoff = 0.2)  #指定q值阈值（可指定 1 以输出全部）
  }

  cluster_up <- data.frame(cluster_up)
  cluster_down <- data.frame(cluster_down)
  write.csv(cluster_up, paste(method, "_up.csv", sep = ""))
  write.csv(cluster_down, paste(method, "_down.csv", sep = ""))

  go_kegg_plot(cluster_up, cluster_down, method, n)
}

