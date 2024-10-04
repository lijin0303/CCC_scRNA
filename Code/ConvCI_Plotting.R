source("~/Rotation/Rafa_Rotation/Analyses/Utilities.R")
sizeSig  <-readRDS("~/Downloads/Latest_Results/SigIPs_3Size.rds")
addSmallLegend <- function(myPlot, pointSize = rel(1.5), textSize = rel(0.55), spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(nrow=8,override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
ErrorBar_IPscore <- function(s,RF="~/Downloads",manu=c(""),out="~/Desktop",
                             w=8,h=7,nr=1){
  data  <- readRDS(glue("{RF}/Latest_Results/ExtracrSeurat_Size{s}.rds"))
  SumIP <- data[,c(1,6:ncol(data))]%>%
    gather(Side_IP,score,-celltype)%>%
    group_by(celltype,Side_IP)%>%
    summarise(cih  = mean(score)+1.96*(sd(score)/sqrt(length(score))),
              cil  = mean(score)-1.96*(sd(score)/sqrt(length(score))),
              mean = mean(score))%>%
    ungroup()
  SumIP%<>%separate(Side_IP,into=c("side","ip"),sep="_")
  KeyIP <- SumIP%>%
    group_by(ip,side)%>%
    arrange(-mean) %>% 
    mutate(rank = 1:n()) %>% 
    filter(rank<4  & grepl("CD8_PD1Tim_",celltype) & side=="ligands")%>%
    pull(ip) %>% unique()
  
  KeyIP <- setdiff(KeyIP,manu)
  pL <- ggplot(data = SumIP%>%
                 filter(side=="ligands" & ip %in% KeyIP),
               aes(celltype, mean, color = celltype)) +
    geom_point(size = rel(2), position=position_dodge(width=0.5)) +
    geom_errorbar(
      aes(ymin = cil, ymax = cih),
      width = 0.1,
      position=position_dodge(width=rel(0.5))) +
    facet_wrap(~ip,nrow=nr)+
    theme_bw()+
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(y="Ligand Scoring")
  pL <- addSmallLegend(pL)
  pR <- ggplot(data = SumIP%>%
                 filter(side=="receptors" & ip %in% KeyIP),
               aes(celltype, mean, color = celltype)) +
    geom_point(size = rel(2), position=position_dodge(width=0.5)) +
    geom_errorbar(
      aes(ymin = cil, ymax = cih),
      width = 0.1,
      position=position_dodge(width=rel(0.5))) +
    facet_wrap(~ip,nrow=nr)+
    theme_bw()+
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(y="Receptor Scoring")
  pR <- addSmallLegend(pR)
  pC <- ggarrange(pL,pR,nrow = 2,common.legend = T)
  ggsave(plot = pC,filename = glue("{out}/Size{s}_IPscore.pdf"),
         width = w,height = h)
}
ErrorBar_IPscore("Small",manu=c("yellow","salmon"))
ErrorBar_IPscore("Medium",w=10)
ErrorBar_IPscore("Large",manu=c("plum1"),w=11)
