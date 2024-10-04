source("~/Rotation/Rafa_Rotation/Analyses/Utilities.R")
NamedVec <- function(D,vec,name){
  nameVec <- D%>%pull(vec)%>%setNames(D%>%pull(name))
  return(nameVec)
}
CCC_Heatmap_Supp <- function(s,RF="~/Downloads",out="~/Desktop"){
  IPl  <- readRDS(glue("{RF}/Latest_Results/{s}_CCIM.rds"))
  pal <- MetPalettes$Redon[[1]]
  AnnotDF <- IPl$PoolType%>%
    mutate(Main = gsub("\\;.*","",Type),
           MajorVote = gsub("\\;.*","",freq))%>%
    mutate(MajorVote = as.numeric(gsub("%","",MajorVote)))%>%
    arrange(Main,-MajorVote)
  AnnotDF%<>%
    mutate(KeyAnnot = case_when(grepl("CD4",Main)~"CD4",
                                grepl("CD8",Main)~"CD8",
                                grepl("Mono",Main)~"Mono|Mac",
                                grepl("DC",Main)~"DC",
                                grepl("Treg",Main)~"Treg",
                                grepl("Native",Main)~"T",
                                grepl("T-Myeloid",Main)~"Myel",
                                TRUE~Main)) 
  AnnotDF%<>% mutate(KeyAnnot = factor(KeyAnnot,levels = unique(AnnotDF$KeyAnnot)))
  N <- length(levels(AnnotDF$KeyAnnot))
  ann <- anno_block(gp = gpar(fill = pal[1:N]),
                    labels = levels(AnnotDF$KeyAnnot),
                    labels_gp = gpar(col = "black", fontsize = rel(7.5)))
  ha = HeatmapAnnotation(foo = ann)
  haR = rowAnnotation(foo =anno_block(gp = gpar(fill = pal[1:N]),
                                      labels = levels(AnnotDF$KeyAnnot),
                                      labels_gp = gpar(col = "black", fontsize = rel(7.5))))
  split = AnnotDF$KeyAnnot
  coul <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  pdf(glue("{out}/Size{s}_CCIM_Heatmap_Supp.pdf"),width=12,height=10)
  p <- ComplexHeatmap::Heatmap(IPl$CCIM[AnnotDF$pool,AnnotDF$pool], 
                          name = "Interaction", 
                          col = coul,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          column_title = glue("For Size {s}: Sender by Receiver"),
                          row_title = NULL,
                          
                          row_names_side = "left",
                          column_names_side = "bottom",
                          column_labels = AnnotDF$Type,
                          row_labels = AnnotDF$Type,
                          row_names_gp = gpar(col = pal[1:N], fontsize = rel(4)),
                          column_names_gp = gpar(col = pal[1:N], fontsize = rel(4)),
                          column_names_rot = 45,
                          
                          top_annotation = ha,
                          right_annotation = haR,
                          column_split =split,
                          row_split =split,
                          border = TRUE)
  draw(p)
  dev.off()
}
CCC_Heatmap_Main <- function(s,RF="~/Downloads",out="~/Desktop"){
  IPl  <- readRDS(glue("{RF}/Latest_Results/{s}_CCIM.rds"))
  pal <- MetPalettes$Redon[[1]]
  AnnotDF <- IPl$PoolType%>%
    mutate(Main = gsub("\\;.*","",Type),
           MajorVote = gsub("\\;.*","",freq))%>%
    mutate(MajorVote = as.numeric(gsub("%","",MajorVote)))%>%
    arrange(Main,-MajorVote)
  AnnotDF%<>%
    mutate(KeyAnnot = case_when(grepl("CD4",Main)~"CD4",
                                grepl("CD8",Main)~"CD8",
                                grepl("Mono",Main)~"Mono|Mac",
                                grepl("DC",Main)~"DC",
                                grepl("Treg",Main)~"Treg",
                                grepl("Native",Main)~"T",
                                grepl("T-Myeloid",Main)~"Myel",
                                TRUE~Main)) 
  PoolAnnot <- NamedVec(AnnotDF,"Main","pool")
  MeanCCC <- IPl$CCIM%>%
    rownames_to_column("senderPool") %>% 
    gather(receiverPool,CCC,-senderPool)%>%
    mutate(SenderT = PoolAnnot[senderPool],
           ReceiverT = PoolAnnot[receiverPool]) %>% 
    group_by(SenderT,ReceiverT) %>% 
    summarise(meanCCC = mean(CCC))%>% 
    spread(ReceiverT,meanCCC)%>% 
    column_to_rownames("SenderT")
  SumAnnot <- AnnotDF%>%
    select(Main,KeyAnnot)%>%
    distinct() %>% 
    arrange(KeyAnnot)
  SumAnnot%<>% mutate(KeyAnnot = factor(KeyAnnot,levels = unique(SumAnnot$KeyAnnot)))
  N <- length(levels(SumAnnot$KeyAnnot))
  ann <- anno_block(gp = gpar(fill = pal[1:N]),
                    labels = levels(SumAnnot$KeyAnnot),
                    labels_gp = gpar(col = "white", fontsize = rel(7.5)))
  ha = HeatmapAnnotation(foo = ann)
  haR = rowAnnotation(foo =anno_block(gp = gpar(fill = pal[1:N]),
                                      labels = levels(SumAnnot$KeyAnnot),
                                      labels_gp = gpar(col = "white", fontsize = rel(7.5))))
  split = SumAnnot$KeyAnnot
 
  breaks <- c(0,5,10,15,20)
  coul <- circlize::colorRamp2(breaks, rev(brewer.pal(n = length(breaks), name ="RdYlBu")))
  
  
  pdf(glue("{out}/Size{s}_CCIM_Heatmap_Main.pdf"),width=12,height=10)
  p <- ComplexHeatmap::Heatmap(as.matrix(MeanCCC[SumAnnot$Main,SumAnnot$Main]), 
                          name = "Interaction", 
                          col = coul,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          column_title = glue("For Size {s}: Sender by Receiver"),
                          row_title = NULL,
                          row_names_side = "left",
                          column_names_side = "bottom",
                          
                          column_labels = SumAnnot$Main,
                          row_labels = SumAnnot$Main,
                          row_names_gp = gpar(col = pal[1:N], fontsize = rel(7.5)),
                          column_names_gp = gpar(col = pal[1:N], fontsize = rel(7.5)),
                          column_names_rot = 45,
                          
                          top_annotation = ha,
                          right_annotation = haR,
                          column_split =split,
                          row_split =split,
                          border = TRUE)
  draw(p)
  dev.off()
}


