source("~/Rotation/Rafa_Rotation/Analyses/Utilities.R")
#===== Key Utility Loading =====
EasyFeature_Ori <- function(data,feature,cols=NULL,scaleColor=1){
  side <- gsub("\\_.*","",feature)
  if(is.null(cols)){
    cols <-  switch(side,
                    "ligands" = c("gray90","#FD8D3C","#BD0026"),
                    "receptors" = c("gray90","#8C96C6","#810F7C"))}
  data%<>%.[,c("cor1","cor2","tnt",feature)]
  data.feature <- data[, feature]
  brewer.gran <- length(cols)*scaleColor
  data$feature <- as.numeric(as.factor(cut(as.numeric(data.feature),
                                           breaks =brewer.gran)))
  TD<- data%>%filter(tnt=="T cells")
  nTD<- data%>%filter(tnt=="Non T cells")
  tfp <- ggplot(TD,aes(x=cor1,y=cor2,color=feature))+
    scale_colour_gradientn(colours = cols)+
    geom_point(data = subset(TD,feature!=max(TD$feature)))+
    geom_point(data = subset(TD,feature==max(TD$feature)))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "right")+
    labs(x="UMAP-1",y="UMAP-2")
  ntfp <- ggplot(nTD,aes(x=cor1,y=cor2,color=feature))+
    scale_colour_gradientn(colours = cols)+
    geom_point(data = subset(nTD,feature!=max(nTD$feature)))+
    geom_point(data = subset(nTD,feature==max(nTD$feature)))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "right")+
    labs(x="UMAP-1",y="UMAP-2")
  fp <- ggarrange(tfp,ntfp,nrow=1,common.legend = T)
  return(fp)
  
}
EasyDim_Ori <- function(data){
  TD<- data%>%filter(tnt=="T cells")
  nTD<- data%>%filter(tnt=="Non T cells")
  tcp <- ggplot(TD,aes(x=cor1,y=cor2,color=celltype))+
    geom_point()+
    geom_label_repel(data = .%>% 
                       group_by(celltype) %>%
                       select(cor1, cor2) %>% 
                       summarize_all(mean),
                     aes(label=celltype),
                     color="black",
                     size = rel(3),
                     max.overlaps=Inf,
                     box.padding = unit(0.1, "lines"))+
    scale_colour_manual(values=met.brewer("Signac",n=length(unique(TD$celltype))))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(x="UMAP-1",y="UMAP-2")
  ntcp <- ggplot(nTD,aes(x=cor1,y=cor2,color=celltype))+
    geom_point()+
    geom_label_repel(data = .%>% 
                       group_by(celltype) %>%
                       select(cor1, cor2) %>% 
                       summarize_all(mean),
                     aes(label=celltype),
                     color="black",
                     size = rel(2.5),
                     max.overlaps=Inf,
                     box.padding = unit(0.1, "lines"))+
    scale_colour_manual(values=met.brewer("Signac",n=length(unique(nTD$celltype))))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(x="UMAP-1",y="UMAP-2")
  cp <- ggarrange(tcp,ntcp,nrow=1)
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
  p1 <- EasyFeature_Ori(data,feature = glue("ligands_{ip}"))
  p2 <- EasyFeature_Ori(data,feature = glue("receptors_{ip}"))
  p3 <- EasyDim_Ori(data)
  pumap <- ggarrange(p1,p3,p2,nrow=3,ncol=1)
  p <- ggarrange(tab,pumap,nrow=2,heights = c(1,3))
  return(p)
}
IP_Size <- function(s,ip,RF="~/Rotation/Rafa_Rotation/Analyses"){
  data  <- readRDS(glue("{RF}/Latest_Results/ExtracrSeurat_Size{s}.rds"))
  cors <- readRDS(glue("{RF}/Latest_Results/Original_Coordinates_umap.rds")) %>% 
    set_colnames(c("cor1","cor2","tnt"))
  data%<>%cbind(cors[rownames(data),])
  li <- tolower(strsplit(s,split ="")[[1]][1])
  ip_details <- readRDS(glue("{RF}/Latest_Results/SigIPs_3Size.rds"))[[li]]
  pc <- LigRecUMAP(data,ip,ip_details)
  return(pc)
}
DynamicIP <- function(lrpair,RF="~/Rotation/Rafa_Rotation/Analyses"){
  ip_details <- readRDS(glue("{RF}/Latest_Results/SigIPs_3Size.rds"))
  sl <- c("small","medium","large")
  allComb <- lapply(seq_along(sl),function(i) ip_details[[i]][,1:2]%>%mutate(tumorSize = sl[i])) %>% 
    bind_rows()
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
              inherit.aes = F,size=rel(3.5))+
    labs(x="",y="")+
    theme_minimal()+
    theme(legend.position = "none")
  
  cccSummary<- lapply(c("Small","Medium","Large"),function(s){
    meanCCC <- readRDS(glue("{RF}/Latest_Results/{s}_CCC.rds"))$PairM[,lrPAIRs$lr_pair]%>%apply(.,2,mean)
    probCCC <- readRDS(glue("{RF}/Latest_Results/{s}_CCC.rds"))$PairM[,lrPAIRs$lr_pair]%>%apply(.,2,function(i) mean(i!=0))
    o <- data.frame(meanCCC=meanCCC,probCCC=probCCC,tumorSize=s)%>%rownames_to_column("pair")
  })%>%bind_rows()%>%
    mutate(tumorSize = factor(tumorSize,levels= c("Small","Medium","Large")))%>%
    mutate(pair = forcats::fct_reorder(pair, meanCCC, .desc = TRUE))
  
  BalloonPlot <- ggplot(cccSummary, aes(x = tumorSize, y = pair)) +
    #geom_point(aes(size = meanCCC,color=probCCC), shape = 21, colour = "black", fill = "cornsilk")+
    geom_point(aes(size = meanCCC,fill=probCCC), shape = 21, colour = "black")+
    scale_size(range = c(0, rel(20)))+
    labs(x="",y="")+
    theme_minimal()+
    scale_fill_gradient2(low = muted("gray"), mid = "lightyellow",high = muted("red"), na.value = NA,breaks=c(0.6,0.8,0.95))
  ZoomPair <- ggarrange(alluvialD,BalloonPlot,nrow=1,widths = c(1.3,1))
  return(ZoomPair)
}
#===== Plotting l=====
ggsave("~/Desktop/Small_ZoomIP_Umap.pdf",IP_Size("Small",ip="darkred"),width = 7.5,height = 15)
ggsave("~/Desktop/Medium_ZoomIP_Umap.pdf",IP_Size("Medium",ip="skyblue"),width = 7.5,height = 15)
ggsave("~/Desktop/Large_ZoomIP_Umap.pdf",IP_Size("Large",ip="orange"),width = 7.5,height = 15)
ggsave("~/Desktop/Large_ZoomIP_Umap2.pdf",IP_Size("Large",ip="violet"),width = 7.5,height = 15)
ggsave("~/Desktop/IP_PairZoom.pdf",DynamicIP("Ccl3=Ccr5",RF="~/Rotation/Rafa_Rotation/Analyses"),width = 12,height = 8.8)
ggsave("~/Desktop/IP_PairZoom.pdf",ZoomPair,width = 12.5,height = 8.8)




