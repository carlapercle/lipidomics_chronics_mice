# Working directory: code

# setwd("your/path/code")  

# example: 
setwd("~/lipidomics/raton/000chronics/code")

rm (list = ls ())

# Exploratory analysis to detect abundance patterns between samples and lipids
# and batch effects anomalous behavior in the data

source("01Functions/function_pcaGenes_2.r")
source("01Functions/plotTreeClust.R")
library(ggplot2)
library(plotly)
library(ggdendro)
library(egg)
library(ggpubr)

# Load data ####
load("00Data/chronic_mice_wo.RData")

samplename <- colnames(df_M_c_wo)

group <- c(substr(colnames(df_M_c_wo), 1, nchar(colnames(df_M_c_wo)) - 1))

sinfo <- as.data.frame(cbind(samplename, group), as.is = T)


palette = c("#264653", "#287271", "#2a9d8f", "#8ab17d", "#e4ba4e", "#f29040", "#f7cab6",
            "#e24d28")  #MetBrewer palette egypt

# Cluster ####
# Euclidean distance

distancia <- dist(t(df_M_c_wo), method = "euclidean")
cluster <- hclust(distancia)

cluster$clase <- sinfo$group

fig_eu_wo <- plotTreeClust(cluster, title = "Clustering, euclidean distance (without outliers)",
                           palette = palette)


# Boxplot ####

muestras <- rep(colnames(df_M_c_wo), each = dim(df_M_c_wo)[1])
valores <- c()
for (columna in colnames(df_M_c_wo)) {
  valores <- c(valores, df_M_c_wo[, columna])
}
rm(columna)

condicion <- rep(group, each = dim(df_M_c_wo)[1])
box_data <- data.frame(muestras, valores, condicion)

box_data2 <- box_data[order(box_data$condicion), ]

f_boxplot <- plot_ly(box_data2, y = ~valores, x = ~muestras, color = ~condicion,
                     colors = palette, type = "box") %>%
  layout(title = "<b> Boxplot by samples </b>", legend = list(title = list(text = "<b> Groups </b>")))

# PCA ####

if(!require('impute')) {
  install.packages('impute')
  library('impute')
}

if(exists(".Random.seed")) rm(.Random.seed)
df.imputed <- impute.knn(as.matrix(df_M_c_wo))


mi.pca <- pcaGenes(df.imputed$data)
mi.pca.df <- as.data.frame(mi.pca$scores)
mi.pca.df$grupo <- sinfo$group
mi.pca.df$var.exp <- round(mi.pca$var.exp * 100, 2)
rownames(mi.pca.df) <- sinfo$samplename

fig_pca <- plot_ly(mi.pca.df, x = ~V1, y = ~V2, z = ~V3, text = rownames(mi.pca.df),
                   color = ~grupo, colors = palette)
fig_pca <- fig_pca %>%
  add_markers(marker = list(size = 4))
fig_pca <- fig_pca %>%
  layout(title = "<b>PCA</b>", legend = list(title = list(text = "<b> Groups </b>")),
         scene = list(xaxis = list(title = paste0("PC3: ", mi.pca.df$var.exp[3], "% variance explained")),
                      zaxis = list(title = paste0("PC1: ", mi.pca.df$var.exp[1], "% variance explained")),
                      yaxis = list(title = paste0("PC2: ", mi.pca.df$var.exp[2], "% variance explained"))))


PC1_PC2 <- ggplot(mi.pca.df, aes(x = V1, y = V2,
                                 colour = grupo, linetype = group, shape = grupo)) +
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  guides(color = guide_legend(title = "Groups"), shape = guide_legend(title = "Groups")) +
  scale_shape_manual(values = c(15, 16,17, 18,15, 16,17, 18)) +
  scale_color_manual(values = palette) +
  geom_point(alpha = 0.8, size = 2.5) +
  stat_ellipse(geom="polygon", aes(fill = group), 
               alpha = 0.2, 
               show.legend = FALSE, 
               level = 0.95) +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent")) +
  xlab(paste("PC1: ", round (mi.pca.df$var.exp[1]),
             "% explained variance", sep = "")) + #Plot PC1 variance
  ylab(paste("PC2: ", round (mi.pca.df$var.exp[2]),
             "% explained variance", sep = ""))


PC3_PC2 <-ggplot(mi.pca.df, aes(x = V3, y = V2,
                                colour = grupo, linetype = group, shape = grupo)) +
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  guides(color = guide_legend(title = "Groups"), shape = guide_legend(title = "Groups")) +
  scale_shape_manual(values = c(15, 16,17, 18,15, 16,17, 18)) +
  scale_color_manual(values = palette) +
  geom_point(alpha = 0.8, size = 2.5) +
  stat_ellipse(geom="polygon", aes(fill = group), 
               alpha = 0.2, 
               show.legend = FALSE, 
               level = 0.95) +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent")) +
  xlab(paste("PC3: ", round (mi.pca.df$var.exp[3]),
             "% explained variance", sep = "")) + #Plot PC1 variance
  ylab(paste("PC2: ", round (mi.pca.df$var.exp[2]),
             "% explained variance", sep = ""))

ggpubr::ggarrange(PC1_PC2,PC3_PC2, ncol= 2, nrow =1, common.legend = T)