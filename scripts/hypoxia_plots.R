library(tidyverse)
library(cowplot)
library(FSA)

setwd("/Users/keithhazleton/Documents/Research/Issa")

pimid_fluor <- read_csv("/Users/keithhazleton/Documents/Research/Issa/pimid_fluor.csv", 
     col_types = cols(mouse_id = col_integer()))

cecum <- filter(pimid_fluor, location == "Cecum")
proximal <- filter(pimid_fluor, location == "Prox_colon")
distal <- filter(pimid_fluor, location == "Dist_colon")

cecum_plot <- ggplot(cecum, aes(diet, fluorescence)) + theme_bw(base_size = 14) + geom_boxplot() +  
    geom_jitter(height = 0, width = 0.2) + xlab('') + guides (x = guide_axis(angle = 45))

prox_plot <- ggplot(proximal, aes(diet, fluorescence)) + theme_bw(base_size = 14) + geom_boxplot() +  
  geom_jitter(height = 0, width = 0.2) + ylab('') + xlab('') + guides (x = guide_axis(angle = 45))

dist_plot <- ggplot(distal, aes(diet, fluorescence)) + theme_bw(base_size = 14) + geom_boxplot() +  
  geom_jitter(height = 0, width = 0.2) + ylab('') + xlab('')+ guides (x = guide_axis(angle = 45))

fluor_fig <- plot_grid(cecum_plot, prox_plot, dist_plot, labels = c("A", "B", "C"), nrow = 1)

save_plot("fluor_fig.png", fluor_fig)

kruskal.test(fluorescence ~ diet, data = cecum)
dunnTest(fluorescence ~ diet, data = cecum)
kruskal.test(fluorescence ~ diet, data = distal)
kruskal.test(fluorescence ~ diet, data = proximal)
dunnTest(fluorescence ~ diet, data = proximal)

ggplot(cecum, aes(diet, fluorescence)) + theme_bw(base_size = 14) + geom_boxplot() +  
  geom_jitter(height = 0, width = 0.2) + xlab('') + guides (x = guide_axis(angle = 45)) + ylab("Hypoxia") + xlab("diet")
  #+  facet_wrap(. ~ batch)

