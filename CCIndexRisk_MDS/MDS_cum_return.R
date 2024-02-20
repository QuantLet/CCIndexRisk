rm(list=ls())
library(qs); library(Rtsne); library(umap);
library(tidyverse); library(lubridate);

################################################################################
############################## Read data #######################################
################################################################################
if (file.exists("Data/portfolio_cum_return.qs")){
  data <- qread("Data/portfolio_cum_return.qs")
}else{
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
  
  ### Calculate cummulative return
  for (i in 2:ncol(data)){
    for (j in 1:nrow(data)){
      if ((i == 2) & (j == 1)){
        data_cum <- as.data.frame(matrix(rep(NA, (ncol(data)-1)*nrow(data)), nrow = nrow(data)))
      }
      data_cum[j,(i-1)] <- prod(1+data[1:j,i], na.rm = T) - 1
    }
  }
  colnames(data_cum) <- colnames(data[-1])
  ### Save cummulative return
  qsave(data_cum,"Data/portfolio_cum_return.qs")
}

################################################################################
############################## Euclid Distance #################################
################################################################################
### Calculate Euclidean distance
d_eucl <- matrix(rep(NA, (ncol(data))^2), nrow = (ncol(data)))
for (i in (1:ncol(data))){
  for (j in (1:ncol(data))){
    d_eucl[i,j] <- dist(t(data[,c(i,j)]), method = "euclidean", p=2)
  }
}
colnames(d_eucl) <- colnames(data)
rownames(d_eucl) <- colnames(data)


################################################################################
### Calculate MDS
set.seed(12345)
MDS_fit <- cmdscale(d_eucl, eig = TRUE, k = 2)
MDS_plot <- data.frame(x = MDS_fit$points[,1],
                       y = MDS_fit$points[,2])


### Divide data by optimalization method
p_type <- rep(NA, nrow(MDS_plot))
p_type[str_detect(rownames(MDS_plot), "ES")] <- "ES"
p_type[str_detect(rownames(MDS_plot), "VAR")] <- "VaR"
p_type[str_detect(rownames(MDS_plot), "PWR")] <- "PWR"
p_type[str_detect(rownames(MDS_plot), "EXP")] <- "EXP"
p_type[str_detect(rownames(MDS_plot), "Index")] <- "Index"
p_type[str_detect(rownames(MDS_plot), "Index_5")] <- "CRIX"
p_type <- factor(as.factor(p_type), levels = c("VaR", "ES", "EXP", "PWR", "Index", "CRIX"))

### Divide data by rebalancing period
p_reb <- rep(NA, nrow(MDS_plot))
p_reb[str_detect(rownames(MDS_plot), "_30")] <- "30"
p_reb[str_detect(rownames(MDS_plot), "_90")] <- "90"
p_reb[str_detect(rownames(MDS_plot), "Index")] <- "Index"
p_reb[str_detect(rownames(MDS_plot), "Index_5")] <- "CRIX"
p_reb <- factor(as.factor(p_reb), levels = c("30", "90", "Index", "CRIX"))

### Divide data by maximum weight
p_max <- rep(NA, nrow(MDS_plot))
p_max[str_detect(rownames(MDS_plot), c("ES_|VAR_|PWR_|EXP_"))] <- "1"
p_max[str_detect(rownames(MDS_plot), c("ES20_|VAR20_|PWR20_|EXP20_"))] <- "0.2"
p_max[str_detect(rownames(MDS_plot), c("ES30_|VAR30_|PWR30_|EXP30_"))] <- "0.3"
p_max[str_detect(rownames(MDS_plot), c("ES50_|VAR50_|PWR50_|EXP50_"))] <- "0.5"
p_max[str_detect(rownames(MDS_plot), "Index")] <- "Index"
p_max[str_detect(rownames(MDS_plot), "Index_5")] <- "CRIX"
p_max <- factor(as.factor(p_max), levels = c("0.2", "0.3", "0.5", "1", "Index", "CRIX"))

### Outlier indices numbers 
ind <- c(rep(NA, 2208), rep(NA, 9), 10, NA, NA, 13, rep(NA, 5), 19, NA, NA, 22)

################################################################################
### Plot MDS by type 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=MDS_plot[which(p_type=="CRIX"),1],y=MDS_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Euclide_distance/cummulative/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by rebalancing period  
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=MDS_plot[which(p_reb=="CRIX"),1],y=MDS_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Euclide_distance/cummulative/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by maximum weight 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=MDS_plot[which(p_max=="CRIX"),1],y=MDS_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Euclide_distance/cummulative/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)



################################################################################
########################### Pearson Correlation ################################
################################################################################
### Calculate Pearson Correlation
d_pear <- cor(data, use="pairwise.complete.obs", method = "pearson")

################################################################################
### Calculate MDS
set.seed(12345)
MDS_fit <- cmdscale(1-d_pear, eig = TRUE, k = 2)
MDS_plot <- data.frame(x = MDS_fit$points[,1],
                       y = MDS_fit$points[,2])

### Outlier indices numbers 
ind <- c(rep(NA, 2208), 1, NA, NA, 4, NA, NA, 7, NA, NA, 10, 11, 12, 13, rep(NA, 7), 21, 22)

################################################################################
### Plot MDS by type 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=MDS_plot[which(p_type=="CRIX"),1],y=MDS_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Pearson_correlation/cummulative/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by rebalancing period  
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=MDS_plot[which(p_reb=="CRIX"),1],y=MDS_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Pearson_correlation/cummulative/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by maximum weight 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=MDS_plot[which(p_max=="CRIX"),1],y=MDS_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Pearson_correlation/cummulative/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)


################################################################################
############################# Cosine Similarity ################################
################################################################################
library(lsa)
### Calculate Cosine similarity
d_cos <- matrix(rep(NA, (ncol(data))^2), nrow = (ncol(data)))
for (i in (1:ncol(data))){
  for (j in (1:ncol(data))){
    M <- as.matrix(na.omit(data[,c(i,j)]))
    d_cos[i,j] <- cosine(x=M[,1], y=M[,2])
  }
}
colnames(d_cos) <- colnames(data)
rownames(d_cos) <- colnames(data)

################################################################################
### Calculate MDS
set.seed(12345)
MDS_fit <- cmdscale(1-d_cos, eig = TRUE, k = 2)
MDS_plot <- data.frame(x = MDS_fit$points[,1],
                       y = MDS_fit$points[,2])

### Outlier indices numbers 
ind <- c(rep(NA, 2208), 1, NA, NA, 4, rep(NA, 5), 10, 11, rep(NA, 9), 21, NA)

################################################################################
### Plot MDS by type 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=MDS_plot[which(p_type=="CRIX"),1],y=MDS_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Cosine_similarity/cummulative/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by rebalancing period  
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=MDS_plot[which(p_reb=="CRIX"),1],y=MDS_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Cosine_similarity/cummulative/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by maximum weight 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=MDS_plot[which(p_max=="CRIX"),1],y=MDS_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/Cosine_similarity/cummulative/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)


################################################################################
################################# L2 distance ##################################
################################################################################
### Calculate L2 distance
d_L2 <- matrix(rep(NA, (ncol(data))^2), nrow = (ncol(data)))
for (i in (1:ncol(data))){
  for (j in (1:ncol(data))){
    d_L2[i,j] <- l2norm(data[,i], data[,j])
  }
}
colnames(d_L2) <- colnames(data)
rownames(d_L2) <- colnames(data)

################################################################################
### Calculate MDS
set.seed(12345)
MDS_fit <- cmdscale(d_L2, eig = TRUE, k = 2)
MDS_plot <- data.frame(x = MDS_fit$points[,1],
                       y = MDS_fit$points[,2])
rownames(MDS_plot) <- colnames(data)

### Outlier indices numbers 
ind <- c(rep(NA, 2208), rep(NA, 21), 22)

################################################################################
### Plot MDS by type 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_type, order = as.numeric(factor(p_type))))+ 
  geom_point(aes(shape=p_type, color=p_type)) +
  geom_point(aes(x=MDS_plot[which(p_type=="CRIX"),1],y=MDS_plot[which(p_type=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("gray", "#000000", "#FF0000",  "#FFA500", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 19, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/L2/cummulative/p_type.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by rebalancing period  
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_reb, order = as.numeric(factor(p_reb))))+ 
  geom_point(aes(shape=p_reb, color=p_reb)) +
  geom_point(aes(x=MDS_plot[which(p_reb=="CRIX"),1],y=MDS_plot[which(p_reb=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/L2/cummulative/p_reb.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)

################################################################################
### Plot MDS by maximum weight 
gg <- ggplot2::ggplot(MDS_plot, aes(x=x, y=y, col=p_max, order = as.numeric(factor(p_max))))+ 
  geom_point(aes(shape=p_max, color=p_max)) +
  geom_point(aes(x=MDS_plot[which(p_max=="CRIX"),1],y=MDS_plot[which(p_max=="CRIX"),2]), color="green", shape=15)+
  geom_text(MDS_plot, mapping=aes(x=x, y=y, label=ind), hjust=0, vjust=0, show.legend = FALSE)+
  theme_classic()+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill='transparent'), legend.background = element_rect(fill='transparent'),
                         plot.background = element_rect(fill='transparent', color=NA))+
  labs(x ="", y = "", color = "", shape = "")+
  xlim(-max(-min(MDS_plot$x), max(MDS_plot$x)), max(-min(MDS_plot$x), max(MDS_plot$x)))+
  ylim(-max(-min(MDS_plot$y), max(MDS_plot$y)), max(-min(MDS_plot$y), max(MDS_plot$y)))+
  scale_color_manual(values = c("#000000", "#FF0000", "#FFA500", "gray", "#0000FF", "#00FF00"))+
  scale_shape_manual(values=c(17, 17, 17, 19, 15, 15))
gg
ggsave("CCIndexRisk_MDS/Output/L2/cummulative/p_max.png", gg, width = 5000, height = 4500, 
       dpi = 500, units = "px", limitsize = F)


