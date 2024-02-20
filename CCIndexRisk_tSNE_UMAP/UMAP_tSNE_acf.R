rm(list=ls())
library(qs); library(Rtsne); library(umap);
library(tidyverse); library(lubridate); library(ggimage);

################################################################################
############################## Read data #######################################
################################################################################

### Portfolio data 90 days a 30 days rebalancing 
data_90 <- qread("Data/portfolio_returns_90d.qs")
colnames(data_90)[2:ncol(data_90)] <- str_c(colnames(data_90)[2:ncol(data_90)], "_90")
data_30 <- qread("Data/portfolio_returns_30d.qs")
colnames(data_30)[2:ncol(data_30)] <- str_c(colnames(data_30)[2:ncol(data_30)], "_30")
data <- cbind(data_90, data_30[,-1])
data$Date <- ymd(data$Date)

### Indices 
temp <- c(list.files("Data/Indices"))

for (k in 1:length(temp)){
  name <- temp[k]
  ind <- read.csv2(paste0("Data/Indices/", name), sep=";", dec=",")
  ind$Date <- dmy(ind$Date) 
  ind <- cbind(ind, "Return" = c(NA, ind[2:nrow(ind),2]/ind[1:(nrow(ind)-1),2]-1))
  colnames(ind)[3] <- paste0("Index_", k)
  
  data <- left_join(data, ind[,c(1,3)], by ="Date")
}

################################################################################
############################# Calculate acf ####################################
################################################################################
for (i in 2:ncol(data)){
  if (i==2){
    data_acf <- acf(data[,i], lag.max = 7, na.action = na.pass)$acf[2:8]
  }else{
    data_acf<- cbind(data_acf, acf(data[,i], lag.max = 7, na.action = na.pass)$acf[2:8])
  }
}

data_acf <- as.data.frame(data_acf)
colnames(data_acf) <- colnames(data[,-1])


################################################################################
################################### tSNE #######################################
################################################################################
### Calculate tSNE
set.seed(12345)
tsne_data <- Rtsne(t(data_acf), dims=2)
tsne_plot <- data.frame(x = tsne_data$Y[,1], y = tsne_data$Y[,2])

### Divide data by optimalization method
p_type <- rep(NA, nrow(tsne_plot))
p_type[str_detect(colnames(data_acf), "ES")] <- "ES"
p_type[str_detect(colnames(data_acf), "VAR")] <- "VaR"
p_type[str_detect(colnames(data_acf), "PWR")] <- "PWR"
p_type[str_detect(colnames(data_acf), "EXP")] <- "EXP"
p_type[str_detect(colnames(data_acf), "Index")] <- "Index"
p_type[str_detect(colnames(data_acf), "Index_5")] <- "CRIX"
p_type <- factor(as.factor(p_type), levels = c("VaR", "ES", "EXP", "PWR", "Index", "CRIX"))

### Divide data by rebalancing period
p_reb <- rep(NA, nrow(tsne_plot))
p_reb[str_detect(colnames(data_acf), "_30")] <- "30"
p_reb[str_detect(colnames(data_acf), "_90")] <- "90"
p_reb[str_detect(colnames(data_acf), "Index")] <- "Index"
p_reb[str_detect(colnames(data_acf), "Index_5")] <- "CRIX"
p_reb <- factor(as.factor(p_reb), levels = c("30", "90", "Index", "CRIX"))

### Divide data by maximum weight
p_max <- rep(NA, nrow(tsne_plot))
p_max[str_detect(colnames(data_acf), c("ES_|VAR_|PWR_|EXP_"))] <- "1"
p_max[str_detect(colnames(data_acf), c("ES20_|VAR20_|PWR20_|EXP20_"))] <- "0.2"
p_max[str_detect(colnames(data_acf), c("ES30_|VAR30_|PWR30_|EXP30_"))] <- "0.3"
p_max[str_detect(colnames(data_acf), c("ES50_|VAR50_|PWR50_|EXP50_"))] <- "0.5"
p_max[str_detect(colnames(data_acf), "Index")] <- "Index"
p_max[str_detect(colnames(data_acf), "Index_5")] <- "CRIX"
p_max <- factor(as.factor(p_max), levels = c("0.2", "0.3", "0.5", "1", "Index", "CRIX"))

### Outlier indices numbers 
ind <- c(rep(NA, 2208), rep(NA, 6), "7", rep(NA, 9), "17", "18", NA, NA, "21", rep(NA,1))

################################################################################
### Plot tSNE by type 
gg <- ggplot2::ggplot(tsne_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=tsne_plot[which(p_type=="CRIX"),1],y=tsne_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(tsne_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(tsne_plot$x), max(tsne_plot$x)), max(-min(tsne_plot$x), max(tsne_plot$x)))+
  ylim(-max(-min(tsne_plot$y), max(tsne_plot$y)), max(-min(tsne_plot$y), max(tsne_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/tSNE/acf/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot tSNE by rebalancing period  
gg <- ggplot2::ggplot(tsne_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=tsne_plot[which(p_reb=="CRIX"),1],y=tsne_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(tsne_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(tsne_plot$x), max(tsne_plot$x)), max(-min(tsne_plot$x), max(tsne_plot$x)))+
  ylim(-max(-min(tsne_plot$y), max(tsne_plot$y)), max(-min(tsne_plot$y), max(tsne_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/tSNE/acf/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot tSNE by maximum weight 
gg <- ggplot2::ggplot(tsne_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=tsne_plot[which(p_max=="CRIX"),1],y=tsne_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(tsne_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(tsne_plot$x), max(tsne_plot$x)), max(-min(tsne_plot$x), max(tsne_plot$x)))+
  ylim(-max(-min(tsne_plot$y), max(tsne_plot$y)), max(-min(tsne_plot$y), max(tsne_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/tSNE/acf/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)


################################################################################
################################### UMAP #######################################
################################################################################
###Calculate UMAP
umap_data <- umap(t(data_acf))
umap_plot <- data.frame(x = as.matrix(umap_data$layout[,1]),
                        y = as.matrix(umap_data$layout[,2]))


### Outlier indices numbers 
ind <- c(rep(NA, 2208), rep(NA, 6), "7", rep(NA, 9), "17", "18", NA, NA, "21", rep(NA,1))

################################################################################
### Plot UMAP by type 
gg <- ggplot2::ggplot(umap_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=umap_plot[which(p_type=="CRIX"),1],y=umap_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(umap_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(umap_plot$x), max(umap_plot$x)), max(-min(umap_plot$x), max(umap_plot$x)))+
  ylim(-max(-min(umap_plot$y), max(umap_plot$y)), max(-min(umap_plot$y), max(umap_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/UMAP/acf/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot UMAP by rebalancing period  
gg <- ggplot2::ggplot(umap_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=umap_plot[which(p_reb=="CRIX"),1],y=umap_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(umap_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(umap_plot$x), max(umap_plot$x)), max(-min(umap_plot$x), max(umap_plot$x)))+
  ylim(-max(-min(umap_plot$y), max(umap_plot$y)), max(-min(umap_plot$y), max(umap_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/UMAP/acf/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot UMAP by maximum weight 
gg <- ggplot2::ggplot(umap_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=umap_plot[which(p_max=="CRIX"),1],y=umap_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(umap_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(umap_plot$x), max(umap_plot$x)), max(-min(umap_plot$x), max(umap_plot$x)))+
  ylim(-max(-min(umap_plot$y), max(umap_plot$y)), max(-min(umap_plot$y), max(umap_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_tSNE_UMAP/Output/UMAP/acf/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

