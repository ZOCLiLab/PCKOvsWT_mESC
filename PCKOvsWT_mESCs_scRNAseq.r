source("./QC_functions.r")

##################################################
##              RUN data pretreatment           ##
##################################################

dat = load_10X_data(workDir %+% "/merged_10X/")
report_dir = workDir 
result = pretreat_cellFiltration(raw_expr_mat=dat, report=TRUE, report_dir=report_dir)
dat_f = pretreat_mitoFiltration(result$expr_mat, report=TRUE, report_dir=report_dir, mito_glist=mtlist, ratio=0.6)
so = CreateSeuratObject(dat_f, project=basename(datasets[di]))
genes=c(nrow(dat),nrow(dat_f))
cells=c(ncol(dat),ncol(dat_f))
risk_genes=result$pc_glist

write.table(colnames(so), workDir %+% "mESC_merged.filtered.all.cellnames.tsv")
write.table(rownames(so), workDir %+% "mESC_merged.filtered.all.genenames.tsv")
saveRDS(so, workDir %+% "/mESC_merged.rds")

rm(dat, dat_f, result, so)
gc()

##################################################
##        RUN SAVER for merged mESC cells       ##
##################################################
require(SAVER)

mESC_merged = loadRDS(workDir %+% "/mESC_merged.rds")

mESC_cells_selected = read.table(workDir %+% "mESC_merged.filtered.all.cellnames.tsv")
mESC_genes_selected = read.table(workDir %+% "mESC_merged.filtered.all.genenames.tsv")

mESC_genes_selected = mESC_genes_selected[-1,1]
mESC_cells_selected = mESC_cells_selected[-1,1]

mESC_expr = mESC_merged@raw.data
mESC_genes_ind = which(rownames(mESC_expr) %in% mESC_genes_selected)
mESC_cells_ind = which(colnames(mESC_expr) %in% mESC_cells_selected)
mESC_expr_f = mESC_expr[mESC_genes_ind,mESC_cells_ind]

rm(mESC_merged)

sections=10
tmp_saver_fileprefix = workDir %+% "/mESC_merged_saver/"
for (i in c(1:10)){
  sgap = round(dim(mESC_expr_f)[1]/10)
  if (i!=10){
    tmp_mESC_merged.saver = saver(as.matrix(mESC_expr_f),
								  ncores=5,
								  pred.genes=c((sgap*(i-1)+1):(sgap*i)),
								  pred.genes.only=TRUE)
								  saveRDS(tmp_mESC_merged.saver,
								  file=paste(tmp_saver_fileprefix,i,".rds"))
    gc()
  }
  else{
    tmp_mESC_merged.saver = saver(as.matrix(mESC_expr_f),
								  ncores=5,
								  pred.genes = c((sgap*(i-1)+1):(dim(mESC_expr_f)[1])),
								  pred.genes.only=TRUE)
								  saveRDS(tmp_mESC_merged.saver,
								  file=paste(tmp_saver_fileprefix,i,".rds"))
    gc()
  }
}
solist = list()
for (i in c(1:10)){
	solist[[i]] = loadRDS(paste(tmp_saver_fileprefix,i,".rids"))
}
saver_so = combine.saver(solist)
rm(solist)
gc()
saver_so = CreateSeuratObject(mESC_merged.saver.all$estimate,project = "10X_PDGFCKOvsWT")
saver_so = NormalizeData(saver_so)
saver_so = ScaleData(saver_so)
# saver_so@imputed = mESC_merged.saver.all$estimate
saver_so@meta.data$group = as.character(unlist(lapply(rownames(saver_so@meta.data),ident_celltype)))
saver_so = FindVariableFeatures(saver_so,
	                             selection.method = "vst",
	                             mean.function=ExpMean,
	                             dispersion.function=LogVMR,
	                             do.plot=TRUE)
saver_so = RunPCA(saver_so,genes.print=TRUE)
saver_so = RunTSNE(saver_so,reduction.type="tsne",dim.embed=3)
saver_so = RunUMAP(saver_so, reduction="pca",n.neighbors = 30)

saver_sofile = paste(tmp_saver_fileprefix,"saver_so.rds",sep="/")
saveRDS(saver_so, file=saver_sofile)
##################################################


##################################################
##       RUN Seurat for individual groups       ##
##################################################

## PCKO group
###################################################
pcko_cell_idx = which(saver_so$group == "PCKO")

pcko_so = CreateSeuratObject(saver_so@assays$RNA@counts[,pcko_cell_idx],
                             project="PCKO",assay="counts",min.cells=0,min.features=0,
                             meta.data=saver_so@meta.data[pcko_cell_idx,])
pcko_so = NormalizeData(pcko_so)
pcko_so = ScaleData(pcko_so)
ccPath = "/home/zhjning/Nutstore/Weisi/mouse_cell_cycle_genes/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds"
mouse.cc.genes = loadRDS(ccPath)
pcko_so = CellCycleScoring(pcko_so, 
                           s.features=mouse.cc.genes$s.genes,
                           g2m.features=mouse.cc.genes$g2m.genes)

pcko_so$CC_diff = pcko_so$S.Score - pcko_so$G2M.Score

pcko_so = ScaleData(object = pcko_so,
                    data.use = pcko_so@assays$counts@scale.data,
                    vars.to.regress = "CC_diff", features = rownames(pcko_so))
pcko_so = FindVariableFeatures(pcko_so)
# pcko_so = RunPCA(pcko_so, features=c(mouse.cc.genes$s.genes, mouse.cc.genes$g2m.genes))
pcko_so = RunPCA(pcko_so)
DimPlot(pcko_so,group.by="Phase")+scale_color_futurama()
pcko_so = FindNeighbors(pcko_so)
pcko_so = RunTSNE(pcko_so,reduction="pca",dims=c(1:20),dim.embed=3,seed.use=999)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.5)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.6)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.8)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=1)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=1.2)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.4)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.3)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.2)
pcko_so = FindClusters(pcko_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.1)
DimPlot(pcko_so,reduction="tsne")
pcko_so = Seurat::RunUMAP(pcko_so,dims=c(1:20),reduction="pca")
p = Seurat::DimPlot(pcko_so,reduction="umap",group.by="counts_snn_res.0.5",label=T,label.size=5)
ggsave(plot=p,filename="~/Desktop/SOT/pcko_dm_umap_counts_snn_res.0.5.pdf",width=6.17,height=4.60,device="pdf")
p = Seurat::DimPlot(pcko_so,reduction="umap",group.by="Phase",label=T,label.size=5)
ggsave(plot=p,filename="~/Desktop/SOT/pcko_dm_umap_Phase.pdf",width=6.17,height=4.60,device="pdf")
p = Seurat::DimPlot(pcko_so,reduction="umap",group.by="RNA_soup_c10n",label=T,label.size=5)
ggsave(plot=p,filename="~/Desktop/SOT/pcko_dm_umap_RNA_soup_c10n.pdf",width=6.17,height=4.60,device="pdf")


Idents(pcko_so) = pcko_so$counts_snn_res.0.4
all_markers_9_cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.3
all_markers_8_cluster = Seurat::FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.2
all_markers_6_cluster = Seurat::FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.1
all_markers_3_cluster = Seurat::FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)

write.table(all_markers_9_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/pcko_9clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_8_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/pcko_8clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_6_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/pcko_6clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_3_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/pcko_3clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)

res="counts_snn_res.0.4"
res="counts_snn_res.0.3"
res="counts_snn_res.0.2"
res="counts_snn_res.0.1"
p = DimPlot(pcko_so,reduction="umap",group.by=res,label=T)
ggsave(plot=p,device="pdf",width=4.09,height=3.25,
       filename=paste0("~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separatred_cluster_plots/pcko_dm_umap_",res,".pdf"))

saveRDS(pcko_so, workDir %+% "/mESC_PCKO.rds")

## WT group
###################################################
wt_cell_idx = which(saver_so$group == "WT")

wt_so = CreateSeuratObject(saver_so@assays$RNA@counts[,wt_cell_idx],
                             project="WT",assay="counts",min.cells=0,min.features=0,
                             meta.data=saver_so@meta.data[wt_cell_idx,])
wt_so = NormalizeData(wt_so)
wt_so = ScaleData(wt_so)
ccPath = "/home/zhjning/Nutstore/Weisi/mouse_cell_cycle_genes/mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds"
mouse.cc.genes = loadRDS(ccPath)
wt_so = CellCycleScoring(wt_so, 
                         s.features=mouse.cc.genes$s.genes,
                         g2m.features=mouse.cc.genes$g2m.genes)

wt_so$CC_diff = wt_so$S.Score - wt_so$G2M.Score

wt_so = ScaleData(object = wt_so,
                  data.use = wt_so@assays$counts@scale.data,
                  vars.to.regress = "CC_diff", features = rownames(wt_so))
wt_so = FindVariableFeatures(wt_so)
# pcko_so = RunPCA(pcko_so, features=c(mouse.cc.genes$s.genes, mouse.cc.genes$g2m.genes))
wt_so = RunPCA(wt_so)
DimPlot(wt_so,group.by="Phase")+scale_color_futurama()
wt_so = FindNeighbors(wt_so,)
wt_so = RunTSNE(wt_so,reduction="pca",dims=c(1:20),dim.embed=3,seed.use=999)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.1)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.2)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.3)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.4)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.5)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.6)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=0.8)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=1)
wt_so = FindClusters(wt_so,random.seed=999,group.singletons=T,n.start=1000,resolution=1.2)
DimPlot(wt_so,reduction="tsne")
wt_so = Seurat::RunUMAP(wt_so,dims=c(1:20),reduction="pca")
p = Seurat::DimPlot(wt_so,reduction="umap",group.by="counts_snn_res.0.5",label=T,label.size=5)

saveRDS(wt_so, workDir %+% "/mESC_WT.rds")

distinct20 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_pal = distinct20[1:13]
color_scheme = data.frame(group=unique(wt_so@meta.data$counts_snn_res.1),color=col_pal)
color_scheme = color_scheme[order(color_scheme$group),]
tsne_3d_loc = wt_so@reductions$tsne@cell.embeddings

Idents(pcko_so) = wt_so$counts_snn_res.1.2
all_markers_14_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.1
all_markers_13_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.8
all_markers_12_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.6
all_markers_11_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.5
all_markers_10_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.4
all_markers_8_2cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.3
all_markers_8_1cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.2
all_markers_6_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = wt_so$counts_snn_res.0.1
all_markers_3_cluster = Seurat::FindAllMarkers(wt_so,slot="counts",only.pos=TRUE)

write.table(all_markers_14_cluster,"~/Desktop/SOT/wt_14clust_markerList.tsv",sep="\t",col.names=T,row.names=T)
write.table(all_markers_13_cluster,"~/Desktop/SOT/wt_13clust_markerList.tsv",sep="\t",col.names=T,row.names=T)
write.table(all_markers_12_cluster,"~/Desktop/SOT/wt_12clust_markerList.tsv",sep="\t",col.names=T,row.names=T)
write.table(all_markers_11_cluster,"~/Desktop/SOT/wt_11_clust_markerList.tsv",sep="\t",col.names=T,row.names=T)
write.table(all_markers_10_cluster,"~/Desktop/SOT/wt_10_clust_markerList.tsv",sep="\t",col.names=T,row.names=T)

write.table(all_markers_8_2cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/wt_8_2clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_8_1cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/wt_8_1clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_6_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/wt_6clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)
write.table(all_markers_3_cluster,"~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separated_cluster_markers/wt_3clust_markerList.tsv",
            sep="\t",col.names=T,row.names=T)

res="counts_snn_res.0.4"

p = DimPlot(wt_so,reduction="umap",group.by=res,label=T)
ggsave(plot=p,device="pdf",width=4.09,height=3.25,
       filename=paste0("~/Desktop/Volumns/Seagate Expansion Drive/Pdgfc_data/archived/Nov2019/group_separatred_cluster_plots/wt_dm_umap_",res,".pdf"))

all_markers_17cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.1
all_markers_16cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.8
all_markers_12cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.6
all_markers_10_1cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)
Idents(pcko_so) = pcko_so$counts_snn_res.0.5
all_markers_10_2cluster = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)

pcko_so@active.ident = pcko_so$pseudo_state
all_markers_pseudo_state = FindAllMarkers(pcko_so,slot="counts",only.pos=TRUE)

saveRDS(pcko_so, workDir %+% "/mESC_PCKO.rds")