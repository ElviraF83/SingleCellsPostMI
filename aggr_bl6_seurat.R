library(Matrix)
require(openxlsx)
library(Seurat)
library(scales)
library(plyr)
library(dplyr)
library(stringr)

setwd("/projects/chenm/scrna_milena/")

set.seed(0)

file_path <- file.path("/projects/furtam/chenm/cellranger_output", "aggr/outs", "filtered_gene_bc_matrices_mex/mm10Zsgreen1")

gene_br_matrix <-  Read10X(data.dir = file_path)
combined <- CreateSeuratObject(raw.data = gene_br_matrix, min.cells=2, min.genes=500)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = combined@data), value = TRUE, ignore.case = TRUE)
percent.mito <- Matrix::colSums(combined@raw.data[mito.genes, ])/Matrix::colSums(combined@raw.data)
combined <- AddMetaData(object = combined, metadata = percent.mito, col.name = "percent.mito")

label<-c("d0", "d1", "d3", "d5", "d7", "d14", "d28") 

foldername <- "aggr_ccregression_1"

ident<-names(combined@ident)
ident[endsWith(ident, "-1")] <- "d0"
ident[endsWith(ident, "-2")] <- "d1"
ident[endsWith(ident, "-3")] <- "d3"
ident[endsWith(ident, "-4")] <- "d5"
ident[endsWith(ident, "-5")] <- "d7"
ident[endsWith(ident, "-6")] <- "d14"
ident[endsWith(ident, "-7")] <- "d28"
combined@meta.data$stim <- factor(ident, levels =c("d0", "d1", "d3", "d5", "d7", "d14", "d28"))


combined <- FilterCells(combined, subset.names = c("nGene", "percent.mito", "nUMI"),
                    low.thresholds = c(200, -Inf, -Inf), high.thresholds=c(5000, 0.1, 15000))

combined <- NormalizeData(combined)
cc.genes <- readLines(con = "/projects/chenm/scrna_milena/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
cc.genes <- str_to_title(cc.genes)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

combined <- CellCycleScoring(object = combined, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
combined@meta.data<-combined@meta.data[, -grep("old.ident",names(combined@meta.data))]
combined <- FindVariableGenes(object = combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 0.5)

combined <- ScaleData(object = combined, vars.to.regress = c("nGene", "nUMI", "percent.mito", "S.Score", "G2M.Score"))

n_dim=10

combined <- RunPCA(object = combined, pc.genes = combined@var.genes)
combined <- RunTSNE(combined,dims.use = 1:n_dim, reduction.use = "pca",do.fast = T)

#plot.genes <- c("Gsn","Dcn", "Dkk3","Apoe","Mt1", "Mt2", "Cthrc1","Acta2","Col1a1","Col1a2",
#                "Saa3","Clu", "C3", "Mfap4")
plot.genes <- c("Gsn","Dcn", "Dkk3","Apoe","Mt1", "Mt2", "Cthrc1","Acta2","Col1a1","Col1a2","Saa3","Clu", "C3", "Mfap4",
                "zsgreen", "Col8a1", "Postn","Cxcl14","Lpl", "Pi16", "Sfrp2")

run_percent<-function(combined, dirname, prefix, res){
  n_cluster = nlevels(combined@ident) 
  my_color_palette <- hue_pal()(n_cluster)
  cell_count = table(combined@ident, combined@meta.data$stim)
  dat2=data.frame(sweep(cell_count,MARGIN=2,FUN="/",STATS=colSums(cell_count)))
  colnames(dat2)=c("cluster", "sample", "percentage")
  dat2 <- ddply(dat2, .(sample),transform, pos = cumsum(percentage) - (0.5 * percentage))
  dat2$sample <- factor(dat2$sample, levels = levels(combined@meta.data$stim))
  p <- ggplot(dat2, aes(x=sample, y=percentage, fill=cluster)) +
    geom_bar(stat='identity',position =position_fill(reverse = TRUE) ,  colour="black") + 
    scale_fill_manual(values=my_color_palette) + scale_y_continuous(expand = c(0,0),
                                                                    limits = c(0,1)) 
  p<-p+ geom_text(data=dat2, aes(x = sample, y = pos, label = paste0(format(round(percentage*100, 1), nsmall =1 ),"%")), size=4)
  
  setEPS()
  postscript(file=paste0(dirname,"/", prefix,'_percentage_res', 10*res, ".eps"), width=6,height=18)
  print(p)
  dev.off()
}

dirname<-file.path("results", foldername)
dir.create(dirname, showWarnings = FALSE)
prefix<-"aggr"

saveRDS(combined, file=paste0(dirname,"/",prefix,"_combined.rds"))

combined1<-combined 
for (res in seq(4,15,1)/10){
    combined <- FindClusters(combined, reduction.type = "pca", 
                              resolution = res)
    p1<-TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "stim") +ggtitle(paste0("resolution = ", res))+theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    p2<-TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)  # + scale_colour_discrete(breaks = ClusterBreaks, labels = ClusterLabels) + labs(x = "t-SNE 1", y = "t-SNE 2")
    p3<-DotPlot(object = combined, genes.plot = plot.genes, plot.legend = TRUE, do.return = T) 
    png(file=paste0(dirname, "/", prefix,"_cluster_res",  10*res, ".png"), width = 800, height = 1600, units = "px")
    print(plot_grid(p1, p2, p3,ncol = 1))
    dev.off()
    run_percent(combined, dirname, prefix, res)
    all_markers <-FindAllMarkers(object = combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
    all_markers<-cbind(all_markers$gene,all_markers[,1:ncol(all_markers)]) 
    all_markers<-split(all_markers, f = all_markers$cluster )
    write.xlsx(all_markers, paste0(dirname,"/", prefix,"_res", 10*res,"_markers.xlsx"), sheetName=paste0("cluster ", 1:length(all_markers)-1))
    combined<-combined1 
}


