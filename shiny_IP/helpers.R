options(stringsAsFactors = F)
if (!require("pacman", quietly = TRUE))
  install.packages("pacman")
if (!require("import", quietly = TRUE))
  install.packages("import")
pacman::p_load(mgsub,dplyr,purrr,readr,openxlsx,igraph,Matrix,
               ggrepel,ggpubr,ggtext,ggplot2,patchwork,ggalluvial,
               RColorBrewer,MetBrewer)
import::from(tidyr,unite,separate,gather,spread,drop_na)
import::from(magrittr,set_colnames,set_rownames,"%<>%")
import::from(tibble,column_to_rownames,rownames_to_column)
import::from(glue,glue)
ResultF <- "05-24"
smallD <- readRDS(glue("Data/{ResultF}/ExtracrSeurat_SizeSmall.rds"))
mediumD <- readRDS(glue("Data/{ResultF}/ExtracrSeurat_SizeMedium.rds"))
largeD <- readRDS(glue("Data/{ResultF}/ExtracrSeurat_SizeLarge.rds"))
sigIP <- readRDS(glue("Data/{ResultF}/SigIPs.rds"))
allComb <- readRDS(glue("Data/{ResultF}/3SizeComb_IP.rds"))

lig_rec_pairs <- unique(allComb$lr_pair)
EasyFeature <- function(data,feature,cols=NULL,scaleColor=1){
  side <- gsub("\\_.*","",feature)
  if(is.null(cols)){
    cols <-  switch(side,
                    "ligands" = c("gray90","#FD8D3C","#BD0026"),
                    "receptors" = c("gray90","#8C96C6","#810F7C"))}
  data%<>%.[,c("UMAP_1","UMAP_2","celltype",feature)]
  data.feature <- data[, feature]
  brewer.gran <- length(cols)*scaleColor
  data$feature <- as.numeric(as.factor(cut(as.numeric(data.feature),
                                           breaks =brewer.gran)))
  fp <- ggplot(data,aes(x=UMAP_1,y=UMAP_2,color=feature))+
    scale_colour_gradientn(colours = cols)+
    geom_point(data = subset(data,feature!=max(data$feature)))+
    geom_point(data = subset(data,feature==max(data$feature)))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(x="UMAP-1",y="UMAP-2")
  return(fp)
  
}
EasyDim <- function(data){
  clusterC <- data %>% 
    group_by(celltype) %>%
    select(UMAP_1, UMAP_2) %>% 
    summarize_all(mean)
  cp <- ggplot(data,aes(x=UMAP_1,y=UMAP_2,color=celltype))+
    geom_point()+
    geom_label_repel(data = clusterC,
                     aes(label=celltype),
                     color="black",
                     size = rel(3.5),
                     max.overlaps=Inf,
                     box.padding = unit(0.5, "lines"))+
    scale_colour_manual(values=met.brewer("Signac",n=length(unique(data$celltype))))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(x="UMAP-1",y="UMAP-2")
  return(cp)
}
LigRecUMAP <- function(data,ip,ip_details){
  IPtab <- ip_details%>%
    filter(name==ip)%>%
    select(lr_pair)%>%
    set_colnames(glue("IP-{ip}"))
  tab <- ggtexttable(IPtab,rows = NULL,
                     theme = ttheme(
                       colnames.style = colnames_style(color = "white", fill = ip,size = rel(15)),
                       tbody.style = tbody_style(color = "black",fill = c("gray90"), 
                                                 size = rel(12))))
  p1 <- EasyFeature(data,feature = glue("ligands_{ip}"))
  p2 <- EasyFeature(data,glue("receptors_{ip}"))
  p3 <- EasyDim(data)
  p <- ggarrange(p1,p2,tab,p3,nrow=2,ncol=2)
  return(p)
}
VlnPlot <- function(data,ip){
  features <- c(glue("ligands_{ip}"),glue("receptors_{ip}"))
  data <- data[,c("celltype",features)]%>%
    set_colnames(c("celltype","Sender","Receiver"))
  inc_lev <- data%>%
    dplyr::group_by(celltype)%>%
    summarise_if(is.numeric, mean, na.rm = TRUE)%>%
    arrange(Sender,Receiver)%>%pull(celltype)
  data$celltype <- factor(
    x = data$celltype,
    levels = inc_lev)
  plist <- lapply(c("Sender","Receiver"),function(f){
    ggplot(data[,c("celltype",f)]%>%
             set_colnames(c("celltype","f")), aes(x=celltype, y=f,fill=celltype)) +
      geom_violin()+
      geom_dotplot(binaxis='y', stackdir='center',binwidth=0.005,dotsize = 0.02)+
      scale_fill_manual(values=met.brewer("Signac",n=30))+
      theme(plot.title = element_text(color='black', hjust = 0.5),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", size = 1, fill = NA),
            panel.grid = element_blank(),
            axis.text.x = element_text(color='black',hjust = 1,angle = 30,
                                       size=rel(0.8)),
            axis.title.x = element_blank(),
            legend.key = element_rect(fill = NA),
            legend.position ="none")+
      labs(title = f,y="module expression")})
  pl <- ggarrange(plotlist = plist,ncol=1)
  return(pl)
}
IP_Size <- function(data,ip,ip_details){
  p1 <- LigRecUMAP(data,ip,ip_details)
  p2 <- VlnPlot(data,ip)
  pc <- ggarrange(p1,p2,ncol=1)
  return(pc)
}
DynamicIP <- function(lrpair){
  allComb%<>%unite("ipSize",name:tumorSize,sep=":",remove=F)
  subIP <- allComb%>%
    filter(lr_pair==lrpair)%>%
    pull(ipSize)
  lrPAIRs <- allComb%>%filter(ipSize %in% subIP)%>%distinct(lr_pair)
  zoomD <- allComb%>%filter(lr_pair %in% lrPAIRs$lr_pair)
  colormap <- unique(zoomD$name)
  names(colormap) = colormap
  zoomD$tumorSize <- factor(zoomD$tumorSize,
                            levels=c("small","medium","large"))
  mapD <- zoomD%>%
    group_by(tumorSize,name)%>%
    summarise(allComb=paste0(lr_pair,collapse = "\n"),
              freq = n())
  mapD$lmax <- ave(mapD$freq, mapD$tumorSize, FUN=function(i){rev(cumsum(rev(i)))})
  mapD$lmin <- ave(mapD$lmax, mapD$tumorSize, FUN=function(i){c(i[-1],0)})
  mapD$yaes <- (mapD$lmax+mapD$lmin)/2
  alluvialD <- ggplot(zoomD,
         aes(x = tumorSize, stratum = name,
             alluvium = lr_pair,fill=name)) +
    scale_fill_manual(values = colormap)+
    geom_flow() +
    geom_stratum(alpha = .35) +
    geom_text(data = mapD,
              aes(y = yaes,
                  x = tumorSize,
                  label = allComb),
              color = "black",
              inherit.aes = F,size=rel(2.5))+
    labs(x="")+
    theme_minimal()+
    theme(legend.position = "none")
  return(alluvialD)
}


