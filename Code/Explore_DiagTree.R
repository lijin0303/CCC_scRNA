setwd("~/Rotation/Rafa_Rotation/Analyses/IP_determination/")
Trees <- readRDS("GeneTree_All.rds")
Dycut <- function(t,h){
  Dymods <- lapply(h,function(i) 
    cutreeDynamic(dendro = t,method="tree",
                  minClusterSize = 5,
                  cutHeight = quantile(t$height, i))%>%
      labels2colors())
  names(Dymods) <- h
  return(Dymods)
}
Hier_Dendro <- function(t,h,DyCut_Out=NULL,s=0.3){
  if(is.null(DyCut_Out)){DyCut_Out = Dycut(t,h)}
  plot(as.dendrogram(t),leaflab = "none",
       main = "Ligand-Receptor Clustering",
       ylab = "Height",ylim=c(min(t$height)-0.02,max(t$height)+0.02))
  colored_bars(bind_cols(DyCut_Out), t, y_shift=s, 
               rowLabels = paste0("Height % = ",h),
               cex.rowLabels=0.5)
  abline(h=h, lty = 5, col="grey36")
  IPList <- lapply(seq_along(h),function(i) 
    data.frame(IP =DyCut_Out[[i]],method = paste("h=",h[i]))%>%
                     dplyr::mutate(pair = t$labels)%>%
                     dplyr::filter(IP!="grey"))%>%
    bind_rows()
  return(list(Dycut = DyCut_Out,IPList=IPList))
}
sL <- Hier_Dendro(Trees$s,h)
mL <- Hier_Dendro(Trees$m,h)
lL <- Hier_Dendro(t=Trees$l,h=seq(0.25,0.6,0.04),s=0.3)
sapply(mL$Dycut,function(s) aricode::ARI(sL$Dycut$`0.48`,s))
colormap <- unlist(sL$Dycut)%>%unique()
names(colormap) = colormap
combIP <- sL$IPList%>%
  filter(method=="h= 0.48")%>%
  mutate(size="small")%>%
  rbind(mL$IPList%>%
          filter(method=="h= 0.48")%>%
          mutate(size="medium"))%>%
  rbind(lL$IPList%>%
         filter(method=="h= 0.49")%>%
         mutate(size="large"))
ggplot(combIP,
       aes(x = size, stratum = IP, alluvium = pair,
           label = IP,fill=IP)) +
  scale_fill_manual(values = colormap)+
  geom_stratum() + 
  geom_flow()+
  #ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 3) +
  labs(x="")+
  theme_minimal() 

ggplot(combIP%>%
         filter(IP %in% c("blue","brown","yellow","turquoise")),
       aes(x = size, stratum = IP, alluvium = pair,
           label = IP,fill=IP)) +
  scale_fill_manual(values = colormap)+
  geom_stratum() + 
  geom_flow()+
  #ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 3) +
  labs(x="")+
  theme_minimal() 

ggplot(combIP%>%
         filter(IP %in% c("brown","yellow","turquoise")),
       aes(x = size, stratum = IP, alluvium = pair,
           label = IP,fill=IP)) +
  scale_fill_manual(values = colormap)+
  geom_stratum() + 
  geom_flow()+
  #ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 3) +
  labs(x="")+
  theme_minimal() 



combIP%>%
  filter(IP %in% c("blue"))%>%
  group_by(pair)%>%
  summarise(n=n())%>%
  filter(n==3)%>%
  arrange(pair)%>%
  select(-n)%>%
  set_rownames(NULL)%>%
  ggtexttable(rows = NULL,theme = ttheme("lBlue"))

combIP%>%
  filter(IP %in% c("turquoise"))%>%
  group_by(pair)%>%
  summarise(n=n())%>%
  filter(n==3)%>%
  arrange(pair)%>%
  select(-n)%>%
  set_rownames(NULL)%>%
  ggtexttable(rows = NULL,theme = ttheme("lGreen"))

combIP%>%
  filter((IP=="brown" &size=="small")|(IP=="brown" &size=="medium"))%>%
  group_by(pair)%>%
  summarise(n=n())%>%
  filter(n==2)%>%
  arrange(pair)%>%
  select(-n)%>%
  set_rownames(NULL)%>%
  ggtexttable(rows = NULL,theme = ttheme("lRed"))

combIP%>%
  filter((IP=="yellow" &size=="small")|(IP=="brown" &size=="medium"))%>%
  group_by(pair)%>%
  summarise(n=n())%>%
  filter(n==2)%>%
  arrange(pair)%>%
  select(-n)%>%
  set_rownames(NULL)%>%
  ggtexttable(rows = NULL,theme = ttheme("lOrange"))

combIP%>%
  filter((IP=="yellow" &size=="large")|(IP=="turquoise" &size=="medium"))%>%
  group_by(pair)%>%
  summarise(n=n())%>%
  filter(n==2)%>%
  arrange(pair)%>%
  select(-n)%>%
  set_rownames(NULL)%>%
  ggtexttable(rows = NULL,theme = ttheme("mOrangeWhite"))


