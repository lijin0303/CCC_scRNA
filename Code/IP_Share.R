rm(list=ls())
setwd("~/Rotation/Rafa_Rotation/Analyses/IPShare_Explore/")
l <- readRDS("IP_details_tumorLarge.rds")
m  <- readRDS("IP_details_tumorMiddle.rds")
s <- readRDS("IP_details_tumorSmall.rds")
abCP <- function(a,b){
  IPs_a <- lapply(a,function(i) strsplit(i,split="\n")[[1]])
  IPs_b <- lapply(b,function(i) strsplit(i,split="\n")[[1]])
  IPab <- lapply(IPs_a,function(s){
    sharedIP <- lapply(IPs_b,intersect,y=s)
    sharedIP <- sharedIP[lengths(sharedIP)>0]
    return(sharedIP)
  })
  IPba <- lapply(IPs_b,function(s){
    sharedIP <- lapply(IPs_a,intersect,y=s)
    sharedIP <- sharedIP[lengths(sharedIP)>(length(s)/2)]
    return(sharedIP)
  })
  IPab <- IPab[lengths(IPab)!=0]
  IPba <- IPba[lengths(IPba)!=0]
  return(list(IPab=IPab,IPba=IPba))
}
test <- abCP(l,m)
## ccl ccr chemokine-related ligand-receptor pairs
## pro: ccl3,4,5,19; ccr5,7,8
## against: ccl21; ccr10
# ICAM1: Endothelial,melanoma
# Reverse for ccr5 and ccr1

