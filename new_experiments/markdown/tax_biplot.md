tax_biplot
================
2024-06-12

``` r
library(rmarkdown)
library(ggpubr)
```

    ## Loading required package: ggplot2

``` r
library(ggplot2)
library(qiime2R)
```

    ## 
    ## Attaching package: 'qiime2R'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     mean_sd

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggh4x)
```

    ## 
    ## Attaching package: 'ggh4x'
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     guide_axis_logticks

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(ggrepel)
```

**Functions**

``` r
prep_my_biplot_files <- function(biplot_fp,
                                 metadata_fp,
                                 tax_fp,
                                 parse_tax = TRUE,
                                 top_tax_n,
                                 dist_pc_axis,
                                 lines_longer_by){
  
  ## reading in file and extracting needed elements for plot
  pcoa_var <- read_qza(biplot_fp)$data$ProportionExplained
  ## for the points on plot
  biplot_points <- read_qza(biplot_fp)$data$Vectors
  ## for the lines w arrows on plot
  biplot_lines <- read_qza(biplot_fp)$data$Species
  
  ## metadata and taxonomy
  meta <- read_tsv(metadata_fp)
  
  ifelse(parse_tax == TRUE, 
         taxonomy <- read_qza(tax_fp)$data %>% 
                      parse_taxonomy() %>%
                      rownames_to_column('asv'),
         taxonomy <- read_qza(tax_fp)$data %>% 
                      rename(asv = Feature.ID)
         )
  ## some data wrangling wooo
  biplot_points <- biplot_points %>% 
                      rename(sampleid = SampleID) %>% 
                      select(sampleid, PC1, PC2, PC3, PC4) %>% 
                      left_join(meta, by = 'sampleid')
  
  biplot_lines <- biplot_lines %>% 
                        mutate(dist_origin1to2 = sqrt(PC1^2 + PC2^2),
                               dist_origin2to3 = sqrt(PC2^2 + PC3^2),
                               dist_origin3to4 = sqrt(PC3^2 + PC4^2)) %>% 
                        top_n(top_tax_n, .data[[dist_pc_axis]]) %>% 
                        mutate(PC1 = (PC1 * 0.3)*lines_longer_by,
                               PC2 = (PC2 * 0.3)*lines_longer_by,
                               PC3 = (PC3 * 0.3)*lines_longer_by,
                               PC4 = (PC4 * 0.3)*lines_longer_by,) %>% 
                        rename(asv = FeatureID) %>% 
                        select(asv, PC1, PC2, PC3, PC4) %>% 
                        left_join(taxonomy, by = 'asv')
  
  ## list of outputs
  my_list <- list(PCvar = pcoa_var,
                  DFPoints = biplot_points,
                  DFLines = biplot_lines,
                  Metadata = meta,
                  Taxonomy = taxonomy)
  return(my_list)
}

## takes the axis variance and makes it into a label for you so you don't
## have to do it every time 
pcoa_ax_lab <- function(unifrac_var, col_name){
  uni_lab <- as.character(round(unifrac_var[col_name] * 100, 2))
  uni_lab <- paste0(col_name, ' - ', uni_lab, '%')
  return(uni_lab)
}

## function for making the actual plot
make_my_biplot <- function(x_axis,
                           y_axis,
                           points_df,
                           points_color,
                           points_shape,
                           color_name,
                           shape_labels,
                           lines_df,
                           lines_label,
                           lines_label_size,
                           x_name,
                           y_name,
                           shape_name,
                           title_content){
     ifelse(class(points_color) == "character", 
         plot <- ggplot() +
                  geom_jitter(data = points_df,
                             aes(x = .data[[x_axis]], y = .data[[y_axis]], 
                                 color = .data[[points_color]],
                                 shape = .data[[points_shape]]),
                             size = 5,
                             alpha = 0.5) +
                 scale_color_viridis(option = 'H', discrete = TRUE,
                                      name = color_name) +
                  scale_shape_discrete(labels = shape_labels) +
                  theme_bw(base_size = 20) +
                  geom_segment(data = lines_df,
                               aes(x = 0, xend = .data[[x_axis]], y = 0, yend = .data[[y_axis]]),
                               arrow = arrow(length = unit(0.3, "cm"))) +
                  geom_text_repel(data = lines_df,
                                  aes(x = .data[[x_axis]], y = .data[[y_axis]], label = .data[[lines_label]]),
                                  color = 'black',
                                  position = 'jitter',
                                  size = lines_label_size,
                                  fontface = 'bold') +
                  labs(x = x_name,
                       y = y_name,
                       shape = shape_name,
                       title = title_content), 
         plot <- ggplot() +
                  geom_jitter(data = points_df,
                             aes(x = .data[[x_axis]], y = .data[[y_axis]],
                                 shape = .data[[points_shape]]),
                             size = 5,
                             alpha = 0.5) +
                 scale_color_viridis(option = 'H', discrete = TRUE,
                                      name = color_name) +
                  scale_shape_discrete(labels = shape_labels) +
                  theme_bw(base_size = 20) +
                  geom_segment(data = lines_df,
                               aes(x = 0, xend = .data[[x_axis]], y = 0, yend = .data[[y_axis]],
                                   color = .data[[lines_label]]),
                               arrow = arrow(length = unit(0.3, "cm")),
                               linewidth = 1.5) +
                  labs(x = x_name,
                       y = y_name,
                       shape = shape_name,
                       title = title_content))
  return(plot)
}
```

**File Paths**

``` r
metadata_FP <- '../data/misc/proc_newExp_d15-d3_metadata.tsv'
tax_FP <- '../data/qiime/taxonomy.qza'
uu_biplot_FP <- '../data/qiime/core_outputs/uu_biplot_matrix.qza'
wu_biplot_FP <- '../data/qiime/core_outputs/wu_biplot_matrix.qza'

vendor_labs <- c('Charles River',
                 'Taconic')
test_biplot <- read_qza(uu_biplot_FP)
```

# **PC1 and PC2**

**Data Wrangling**

``` r
## unweighted unifrac
uu_files <- prep_my_biplot_files(biplot_fp = uu_biplot_FP,
                                 metadata_fp = metadata_FP,
                                 tax_fp = tax_FP,
                                 parse_tax = TRUE,
                                 top_tax_n = 10,
                                 dist_pc_axis = 'dist_origin1to2',
                                 lines_longer_by = 4)
```

    ## Rows: 98 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): sampleid, mouse_sex, diet, vendor, experiment, experiment_set, samp...
    ## dbl (7): study, experiment_day, day_post_inf, high_fat, high_fiber, purified...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
uu_biplot_points <- uu_files$DFPoints %>% 
                      mutate(day_post_inf = as.factor(day_post_inf))

uu_biplot_lines <- uu_files$DFLines
uu_pcoa_var <- uu_files$PCvar

## making x and y axis labels
uu_x_lab <- pcoa_ax_lab(unifrac_var = uu_pcoa_var,
                        col_name = 'PC1')
uu_y_lab <- pcoa_ax_lab(unifrac_var = uu_pcoa_var,
                        col_name = 'PC2')
```

``` r
## weighted unifrac 
wu_files <- prep_my_biplot_files(biplot_fp = wu_biplot_FP,
                                 metadata_fp = metadata_FP,
                                 tax_fp = tax_FP,
                                 parse_tax = TRUE,
                                 top_tax_n = 10,
                                 dist_pc_axis = 'dist_origin1to2',
                                 lines_longer_by = 9)
```

    ## Rows: 98 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): sampleid, mouse_sex, diet, vendor, experiment, experiment_set, samp...
    ## dbl (7): study, experiment_day, day_post_inf, high_fat, high_fiber, purified...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
wu_biplot_points <- wu_files$DFPoints %>% 
                      mutate(day_post_inf = as.factor(day_post_inf))

wu_biplot_lines <- wu_files$DFLines
wu_pcoa_var <- wu_files$PCvar

## making x and y axis labels
wu_x_lab <- pcoa_ax_lab(unifrac_var = wu_pcoa_var,
                        col_name = 'PC1')
wu_y_lab <- pcoa_ax_lab(unifrac_var = wu_pcoa_var,
                        col_name = 'PC2')
```

**Plots**

``` r
## unweighted unifrac 
uu_biplot <- make_my_biplot(x_axis = 'PC1',
                            y_axis = 'PC2',
                            points_df = uu_biplot_points,
                            points_color = 'day_post_inf',
                            points_shape = 'vendor',
                            color_name = 'Days Relative\nto Infection',
                            shape_labels = vendor_labs,
                            lines_df = uu_biplot_lines,
                            lines_label = 'Family',
                            lines_label_size = 5,
                            x_name = uu_x_lab,
                            y_name = uu_y_lab,
                            shape_name = 'Vendor',
                            title_content = 'Unweighted UniFrac Biplot')

## weighted unifrac
wu_biplot <- make_my_biplot(x_axis = 'PC1',
                            y_axis = 'PC2',
                            points_df = wu_biplot_points,
                            points_color = 'day_post_inf',
                            points_shape = 'vendor',
                            color_name = 'Days Relative\nto Infection',
                            shape_labels = vendor_labs,
                            lines_df = wu_biplot_lines,
                            lines_label = 'Family',
                            lines_label_size = 5,
                            x_name = wu_x_lab,
                            y_name = wu_y_lab,
                            shape_name = 'Vendor',
                            title_content = 'Weighted UniFrac Biplot')

uu_biplot
```

![](tax_biplot_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
wu_biplot
```

![](tax_biplot_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

# **PC2 and PC3**

**Data Wrangling**

``` r
## unweighted unifrac
uu_files2_3 <- prep_my_biplot_files(biplot_fp = uu_biplot_FP,
                                 metadata_fp = metadata_FP,
                                 tax_fp = tax_FP,
                                 top_tax_n = 10,
                                 dist_pc_axis = 'dist_origin2to3',
                                 lines_longer_by = 4)
```

    ## Rows: 98 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): sampleid, mouse_sex, diet, vendor, experiment, experiment_set, samp...
    ## dbl (7): study, experiment_day, day_post_inf, high_fat, high_fiber, purified...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
uu_biplot_points2_3 <- uu_files2_3$DFPoints %>% 
                      mutate(day_post_inf = as.factor(day_post_inf))

uu_biplot_lines2_3 <- uu_files2_3$DFLines
uu_pcoa_var2_3 <- uu_files2_3$PCvar

## making x and y axis labels
uu_x_lab2_3 <- pcoa_ax_lab(unifrac_var = uu_pcoa_var2_3,
                           col_name = 'PC2')
uu_y_lab2_3 <- pcoa_ax_lab(unifrac_var = uu_pcoa_var2_3,
                           col_name = 'PC3')
```

``` r
## weighted unifrac 
wu_files2_3 <- prep_my_biplot_files(biplot_fp = wu_biplot_FP,
                                   metadata_fp = metadata_FP,
                                   tax_fp = tax_FP,
                                   top_tax_n = 10,
                                   dist_pc_axis = 'dist_origin2to3',
                                   lines_longer_by = 9)
```

    ## Rows: 98 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): sampleid, mouse_sex, diet, vendor, experiment, experiment_set, samp...
    ## dbl (7): study, experiment_day, day_post_inf, high_fat, high_fiber, purified...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
wu_biplot_points2_3 <- wu_files2_3$DFPoints %>% 
                        mutate(day_post_inf = as.factor(day_post_inf))

wu_biplot_lines2_3 <- wu_files2_3$DFLines
wu_pcoa_var2_3 <- wu_files2_3$PCvar

## making x and y axis labels
wu_x_lab2_3 <- pcoa_ax_lab(unifrac_var = wu_pcoa_var2_3,
                           col_name = 'PC2')
wu_y_lab2_3 <- pcoa_ax_lab(unifrac_var = wu_pcoa_var2_3,
                           col_name = 'PC3')
```

**Plots**

``` r
## unweighted unifrac 
uu_biplot2_3 <- make_my_biplot(x_axis = 'PC2',
                              y_axis = 'PC3',
                              points_df = uu_biplot_points2_3,
                              points_color = 'day_post_inf',
                              points_shape = 'vendor',
                              color_name = 'Days Relative\nto Infection',
                              shape_labels = vendor_labs,
                              lines_df = uu_biplot_lines2_3,
                              lines_label = 'Family',
                              lines_label_size = 5,
                              x_name = uu_x_lab2_3,
                              y_name = uu_y_lab2_3,
                              shape_name = 'Vendor',
                              title_content = 'Unweighted UniFrac Biplot (PC2-PC3)')

## weighted unifrac
wu_biplot2_3 <- make_my_biplot(x_axis = 'PC2',
                              y_axis = 'PC3',
                              points_df = wu_biplot_points2_3,
                              points_color = 'day_post_inf',
                              points_shape = 'vendor',
                              color_name = 'Days Relative\nto Infection',
                              shape_labels = vendor_labs,
                              lines_df = wu_biplot_lines2_3,
                              lines_label = 'Family',
                              lines_label_size = 5,
                              x_name = wu_x_lab2_3,
                              y_name = wu_y_lab2_3,
                              shape_name = 'Vendor',
                              title_content = 'Weighted UniFrac Biplot (PC2-PC3)')

uu_biplot2_3
```

![](tax_biplot_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
wu_biplot2_3
```

![](tax_biplot_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

# **Saving my Outputs**

``` r
ggsave('../plots/uu_biplot.pdf',
       plot = uu_biplot,
       width = 15,
       height = 10)
ggsave('../plots/wu_biplot.pdf',
       plot = wu_biplot,
       width = 15,
       height = 10)
```
