```R
suppressPackageStartupMessages({

    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(repr)
    library(pheatmap)
    library(irr)
    library(lmerTest)
    library(gridExtra)
    library(pROC)
    library(countsplit)
    library(caret)
    
})

set.seed(42)
```

### The original data from paper
Reading the count file for microarray data from [This paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2001-2-10-research0042) containing 24 samples from 4 mice with two differnt starins and 6 brain regions.




```R
cnt <- fread("./Alldatanomef.txt")

gene_ids <- cnt$clone

cnt <- cnt[, 2:25]
cnt <- cnt[, lapply(.SD, as.numeric)]
cnt <- as.matrix(cnt)

rownames(cnt) <- gene_ids
```


```R
meta <- data.table(region = sub("(\\d+)([A-Za-z]+)\\(.*\\.txt", "\\2", colnames(cnt)),
                    strain = sub("(\\d+)([A-Za-z]+)\\(.*\\.txt", "\\1", colnames(cnt)))
meta[, index := .I]
meta$region <- sub("^B(.*$)", "\\1", meta$region)


cell_line_to_keep <- c("Amygdala", "Cerebellum", "Cortex", "EntorhinalCortex", "Hippocampus", "Midbrain")
meta <- meta[region %in% cell_line_to_keep, ]
```

Let's see the original raw data to see if it is normalized or needs normalization.


```R
options(repr.plot.width = 7)
boxplot(cnt)
```


    
![png](output_5_0.png)
    



```R
# measuering the varince per each column of raw matrix
var_asses <- apply(cnt, 2, var) %>% as.data.table(keep.rownames = T)
setnames(var_asses, old = c("rn", "."), new = c("og_colnames", "Var"))

# temp meta for variance assessmnet
meta_assess <- meta
meta_assess$og_colnames <- colnames(cnt) 
var_asses <- merge(var_asses, meta_assess, by = "og_colnames")

# converting the region and strain to factor to be used in ggplot2
var_asses$region <- as.factor(var_asses$region)
var_asses$strain <- as.factor(var_asses$strain)

# adding the ID for each cluster
var_asses$mouse <- as.factor(c(rep(c("mc_1", "mc_2"), 6), rep(c("mc_3", "mc_4"), 6)))
meta$mouse <- as.factor(c(rep(c("mc_1", "mc_2"), 6), rep(c("mc_3", "mc_4"), 6)))
```


```R
options(repr.plot.height = 7, repr.plot.width = 10)
ggplot(data = var_asses) +
      geom_boxplot(aes(y = Var, x = region)) +  
      geom_point(aes(y = Var, x = region, color = region, shape = strain), size = 3) + 
      labs(title = "Variance Between Regions", subtitle = "Divided by Strains", x = "", y = "") + 
      facet_wrap(~strain) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

```


    
![png](output_7_0.png)
    



```R
options(repr.plot.height = 7, repr.plot.width = 11)
p1 <- ggplot(data = var_asses) +
      geom_boxplot(aes(y = Var, x = mouse)) +
      geom_point(aes(y = Var, x = mouse, color = region, shape = mouse), size = 3) + 
      labs(title = "Variance Between Clusters", subtitle = "Divided by Strains", x = "", y = "") + 
      theme_bw()

p2 <- ggplot(data = var_asses) +
      geom_boxplot(aes(y = Var, x = strain)) +
      geom_point(aes(y = Var, x = strain, color = region, shape = mouse), size = 3) + 
      labs(title = "Variance Between Clusters", subtitle = "Divided by Strains", x = "", y = "") + 
      theme_bw()


grid.arrange(p1, p2, ncol = 2)
```


    
![png](output_8_0.png)
    


### Making the template for feature selection


```R
test_fun <- function(new_col_name, criteria, target_col_name, metadata) {
  
  # Initiating the column with zeros
  metadata[, (new_col_name) := 0]
  
  # Finding indexes and updating the new column
  sub_indexes <- which(metadata[[target_col_name]] == criteria)
  
  metadata[sub_indexes, (new_col_name) := 1]
  
  return(metadata)
}
```


```R
test_fun(new_col_name = "starin_temp", criteria = "B6", target_col_name = "strain", metadata = meta)


regions <- unique(meta$region)

# Loop over each region to create binary columns
for (region in regions) {
  column_name <- paste0(region, "_temp")
  meta <- test_fun(column_name, region, "region", meta)
}

```


```R
meta[c(2,5,6,9,11), ]
```


<table class="dataframe">
<caption>A data.table: 5 × 11</caption>
<thead>
	<tr><th scope=col>region</th><th scope=col>strain</th><th scope=col>index</th><th scope=col>mouse</th><th scope=col>starin_temp</th><th scope=col>Amygdala_temp</th><th scope=col>Cerebellum_temp</th><th scope=col>Cortex_temp</th><th scope=col>EntorhinalCortex_temp</th><th scope=col>Hippocampus_temp</th><th scope=col>Midbrain_temp</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Amygdala   </td><td>129</td><td> 2</td><td>mc_2</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><td>Cortex     </td><td>129</td><td> 5</td><td>mc_1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><td>Cortex     </td><td>129</td><td> 6</td><td>mc_2</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><td>Hippocampus</td><td>129</td><td> 9</td><td>mc_1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><td>Midbrain   </td><td>129</td><td>11</td><td>mc_1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td></tr>
</tbody>
</table>



#### Strain-specific expression
***
At a p-value cut-off of 10-5, we found about 65 probe sets that show strain-specific variation in expression across all areas, compared to 24 genes showing overall differences between strains identified by Sandberg et al. [1]. The genes near the top of our ranking include many of those chosen by Sandberg [12], and all 24 genes identified by Sandberg et al. have p-values less than 10-3 by either ANOVA, template match, or both. Two examples of strain-specific genes not identified in the previous work are shown in more detail Figure 4a. The first is Sparc/osteonectin (testican), which was detected at lower levels in the C57B16 strain. The second is phosphatase ACP1/ACP2, with the opposite expression pattern. We rank ACP1/ACP2 and Sparc/osteonectin as having the 9th and 23rd strongest strain differences, respectively (both ACP1/ACP2 and Sparc/osteonectin are also shown in Figure 3). For comparison, Figure 4b shows β-globin, which was identified in [1] and was ranked 67th by template match. Overall, we rank 46 previously unrecognized genes as having stronger strain differences than β-globin.


```R
fit_strain_temp <- lapply(1:nrow(cnt), function(gene){fit_summary <- lm(cnt[gene,] ~ meta$starin_temp)})
```


```R
summary(fit_strain_temp[[2]])
```


    
    Call:
    lm(formula = cnt[gene, ] ~ meta$starin_temp)
    
    Residuals:
         Min       1Q   Median       3Q      Max 
    -111.417  -29.167    1.917   33.583   85.583 
    
    Coefficients:
                     Estimate Std. Error t value Pr(>|t|)    
    (Intercept)       172.750     13.947  12.386 2.16e-11 ***
    meta$starin_temp   -7.333     19.724  -0.372    0.714    
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    
    Residual standard error: 48.31 on 22 degrees of freedom
    Multiple R-squared:  0.006244,	Adjusted R-squared:  -0.03893 
    F-statistic: 0.1382 on 1 and 22 DF,  p-value: 0.7136




```R
gene_2_dt <- data.table(gene = as.numeric(cnt[2,]), region = factor(meta$region), starin = factor(meta$strain), mouse = meta$mouse)
```


```R
options(repr.plot.height = 7, repr.plot.width = 14)

p1 <- ggplot(gene_2_dt) +
geom_point(aes(x = mouse, y = gene, color = region, shape = starin), size = 3) +
theme_bw() +
geom_abline(slope = -7.333, intercept = 172.750, color = "red") +
labs(x = "Regions", y = "Count", title = "Gene 2 conventional linear model fit", subtitle = "Divided by mice samples")


p2 <- ggplot(gene_2_dt) +
geom_point(aes(x = starin, y = gene, color = region, shape = mouse), size = 3) +
theme_bw() +
geom_abline(slope = -7.333, intercept = 172.750, color = "red") +
labs(x = "Regions", y = "Count", title = "Gene 2 conventional linear model fit", subtitle = "Divided by mice starin (used in the model)")


grid.arrange(p2, p1, ncol = 2)
```


    
![png](output_17_0.png)
    



```R
mixed_effect_model <- lmerTest::lmer(cnt[2,] ~ meta$starin_temp + (meta$starin_temp | meta$mouse))
#summary(mixed_effect_model)
ranef(mixed_effect_model)
```

    boundary (singular) fit: see help('isSingular')
    



    $`meta$mouse`
           (Intercept) meta$starin_temp
    mc_1  5.967905e+00    -5.967858e+00
    mc_2 -5.967905e+00     5.967858e+00
    mc_3 -2.776808e-05     2.774694e-05
    mc_4  2.776808e-05    -2.774694e-05
    
    with conditional variances for “meta$mouse” 



```R
p1 <- ggplot(gene_2_dt) +
geom_point(aes(x = starin, y = gene, color = region, shape = mouse), size = 3) +
theme_bw() +
geom_abline(slope = -7.333 , intercept = 172.750 - 5.967858e+00, color = "blue", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 + 5.967858e+00, color = "orange", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 + 2.774694e-05, color = "green", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 - 2.774694e-05, color = "gray", linetype = 5) + 
labs(x = "Regions", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Divided by mice starin (used in the model)")



p2 <- ggplot(gene_2_dt) +
geom_point(aes(x = mouse, y = gene, color = region, shape = starin), size = 3) +
theme_bw() +
geom_abline(slope = -7.333 , intercept = 172.750 - 5.967858e+00, color = "blue", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 + 5.967858e+00, color = "orange", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 + 2.774694e-05, color = "green", linetype = 5) + 
geom_abline(slope = -7.333, intercept = 172.750 - 2.774694e-05, color = "gray", linetype = 5) + 
labs(x = "Regions", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Divided by mice samples")



grid.arrange(p1, p2, ncol = 2)
```


    
![png](output_19_0.png)
    



```R
options(repr.plot.height = 10, repr.plot.width = 11)


p1 <- ggplot(gene_2_dt[mouse == "mc_1", ]) +
geom_point(aes(x = region, y = gene, color = region, shape = starin), size = 3) +
theme_bw() + 
geom_abline(slope = -7.333 , intercept = 172.750 + 5.967858e+00, color = "#10c470", linetype = 5)  + 
geom_abline(slope = -7.333 , intercept = 172.750 , color = "red", linetype = 1) + 
labs(x = "", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Mouse 1") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))



p2 <- ggplot(gene_2_dt[mouse == "mc_2", ]) +
geom_point(aes(x = region, y = gene, color = region, shape = starin), size = 3) +
theme_bw() + 
geom_abline(slope = -7.333 , intercept = 172.750 - 5.967905e+00, color = "#10c470", linetype = 5)  + 
geom_abline(slope = -7.333 , intercept = 172.750 , color = "red", linetype = 1) + 
labs(x = "", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Mouse 2") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))



p3 <- ggplot(gene_2_dt[mouse == "mc_3", ]) +
geom_point(aes(x = region, y = gene, color = region, shape = starin), size = 3) +
theme_bw() + 
geom_abline(slope = -7.333 , intercept = 172.750 - 2.776808e-05, color = "#10c470", linetype = 5)  + 
geom_abline(slope = -7.333 , intercept = 172.750 , color = "red", linetype = 1) + 
labs(x = "", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Mouse 3") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))


p4 <- ggplot(gene_2_dt[mouse == "mc_4", ]) +
geom_point(aes(x = region, y = gene, color = region, shape = starin), size = 3) +
theme_bw() + 
geom_abline(slope = -7.333 , intercept = 172.750 + 2.776808e-05, color = "#10c470", linetype = 5)  + 
geom_abline(slope = -7.333 , intercept = 172.750 , color = "red", linetype = 1) + 
labs(x = "", y = "Count", title = "Gene 2 mixed effect linear model fit", subtitle = "Mouse 4") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
```


    
![png](output_20_0.png)
    



```R
p_values_strain_temp <- sapply(fit_strain_temp, function(fit) summary(fit)$coefficients[2, 4])
```


```R
r_squ_strain_temp <- sapply(fit_strain_temp, function(fit) summary(fit)$r.squared)
```


```R
assess_df <- data.table(p_val = p_values_strain_temp, 
                       r_squ = r_squ_strain_temp,
                       features = as.factor(rownames(cnt)))
```


```R
temp_res <- data.table(genes = rownames(cnt), p_val = p_values_strain_temp)
temp_res <- temp_res[p_val < 1e-5, ]
```


```R
selected_genes <- temp_res[order(p_val), ] 
```


```R
index_cnt <- which(rownames(cnt) %in% selected_genes$genes)
cnt_temp <- cnt[index_cnt, ]
```


```R
options(repr.plot.height = 6, repr.plot.width = 12)

p1 <- ggplot(data =NULL) +
        geom_histogram(aes(x = r_squ_strain_temp), bins = 10, color = "black", fill = "#788ab3") +
        theme_bw() +
        labs(x = "R^2", y = "", title = "Coefficient of Determination", subtitle = "Conventional Method - Entire Dataset")


p2 <- ggplot(data = NULL) +
        geom_histogram(aes(x = p_values_strain_temp), bins = 10, color = "black", fill = "#8b78b3") +
        theme_bw() +
        labs(x = "P values", y = "", title = "Distribution of P Values", subtitle = "Conventional Method - Entire Dataset")


grid.arrange(p1, p2, ncol = 2)
```


    
![png](output_27_0.png)
    



```R
p1 <- ggplot(data = assess_df[index_cnt, ]) +
geom_histogram(aes(x = r_squ), bins = 10, color = "black", fill = "#788ab3") +
theme_bw() +
labs(x = "R^2", y = "", title = "Coefficient of Determination", subtitle = "Conventional Method - Significant Hits") +
xlim(c(0,1))



p2 <- ggplot(data = assess_df[index_cnt, ]) +
geom_histogram(aes(x = p_val), bins = 10, color = "black", fill = "#8b78b3") +
theme_bw() +
labs(x = "P values", y = "", title = "Distribution of P Values", subtitle = "Conventional Method - Significant Hits")
grid.arrange(p1, p2, ncol = 2)
```

    Warning message:
    “[1m[22mRemoved 2 rows containing missing values or values outside the scale range
    (`geom_bar()`).”



    
![png](output_28_1.png)
    



```R
options(repr.plot.height = 7, repr.plot.width = 15)
pheatmap(t(cnt_temp), scale = "column", main = "Strain Dependent")
```


    
![png](output_29_0.png)
    


The dataset of Sandberg et al. [1] consists of duplicate analysis of six brain regions (**amygdala, cerebellum, cortex, entorhinal cortex, hippocampus and midbrain**) in two strains of mice (129SvEv and C57BL/6), for 13,067 genes and ESTs. We performed two-factor ANOVA and feature selection to look for strain- and/or region-specific variation in gene expression in this data. Our feature-selection strategy, which we call 'template matching', is depicted schematically in Figure 1.

## Mixed effect models
**Hint**:This step might took a while, somewhere between 1-30 minutes depending on your compute power.


```R
suppressWarnings({
    suppressMessages({

        # the random effect is so small, that's why we are getting the `boundary (singular) fit: see help('isSingular')` message
        fit_starin_temp_MEM <- lapply(1:nrow(cnt), function(gene){fit_summary <- lmerTest::lmer(cnt[gene,] ~ meta$starin_temp + (meta$starin_temp | meta$mouse))}) 
})
    })
```


```R
p_values_fit_starin_temp_MEM <- sapply(fit_starin_temp_MEM, function(fit) summary(fit)$coefficients[2, 5])
```


```R
options(repr.plot.height = 7, repr.plot.width = 12)

p1 <- ggplot(data =NULL) +
        geom_histogram(aes(x = p_values_strain_temp), bins = 10, color = "black", fill = "#b38378") +
        theme_bw() +
        labs(x = "P values", y = "", title = "Distribution of P Values", subtitle = "Conventional Method - Entire Dataset") +
        ylim(c(0, 2000))



p2 <- ggplot(data = NULL) +
geom_histogram(aes(x = p_values_fit_starin_temp_MEM), bins = 10, color = "black", fill = "#78b38d") +
theme_bw() +
labs(x = "P values", y = "", title = "Distribution of P Values", subtitle = "Mixed effect method - Entire Dataset") + 
ylim(c(0, 2000))


grid.arrange(p1, p2, ncol = 2)
```


    
![png](output_34_0.png)
    



```R
temp_res <- data.table(genes = rownames(cnt), p_val = p_values_fit_starin_temp_MEM)
temp_res <- temp_res[p_val < 1e-5, ]

selected_genes <- temp_res[order(p_val), ] 

index_cnt <- which(rownames(cnt) %in% selected_genes$genes)
cnt_temp <- cnt[index_cnt, ]

options(repr.plot.height = 7, repr.plot.width = 15)
pheatmap(t(cnt_temp), scale = "column", main = "Strain Dependent")
```


    
![png](output_35_0.png)
    



```R
sessionInfo()
```


    R version 4.3.1 (2023-06-16)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: CentOS Linux 7 (Core)
    
    Matrix products: default
    BLAS/LAPACK: FlexiBLAS IMKL;  LAPACK version 3.11.0
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: America/New_York
    tzcode source: system (glibc)
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] caret_6.0-94       lattice_0.21-8     countsplit_4.0.0   pROC_1.18.5       
     [5] performance_0.12.4 gridExtra_2.3      lmerTest_3.1-3     lme4_1.1-35.3     
     [9] Matrix_1.6-5       irr_0.84.1         lpSolve_5.6.21     pheatmap_1.0.12   
    [13] repr_1.1.7         dplyr_1.1.4        ggplot2_3.5.1      data.table_1.15.4 
    
    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1     timeDate_4041.110    IRdisplay_1.1.0.9000
     [4] farver_2.1.2         fastmap_1.2.0        digest_0.6.36       
     [7] rpart_4.1.19         timechange_0.3.0     lifecycle_1.0.4     
    [10] Cairo_1.6-2          survival_3.5-5       magrittr_2.0.3      
    [13] compiler_4.3.1       rlang_1.1.4          tools_4.3.1         
    [16] utf8_1.2.4           labeling_0.4.3       plyr_1.8.9          
    [19] RColorBrewer_1.1-3   pbdZMQ_0.3-11        withr_3.0.0         
    [22] purrr_1.0.2          numDeriv_2016.8-1.1  stats4_4.3.1        
    [25] nnet_7.3-19          grid_4.3.1           fansi_1.0.6         
    [28] colorspace_2.1-1     future_1.33.2        globals_0.16.3      
    [31] scales_1.3.0         iterators_1.0.14     MASS_7.3-60         
    [34] insight_0.20.5       cli_3.6.3            crayon_1.5.3        
    [37] generics_0.1.3       future.apply_1.11.2  reshape2_1.4.4      
    [40] minqa_1.2.6          stringr_1.5.1        splines_4.3.1       
    [43] parallel_4.3.1       base64enc_0.1-4      vctrs_0.6.5         
    [46] hardhat_1.4.0        boot_1.3-28.1        jsonlite_1.8.8      
    [49] listenv_0.9.1        foreach_1.5.2        gower_1.0.1         
    [52] recipes_1.1.0        glue_1.7.0           parallelly_1.37.1   
    [55] nloptr_2.0.3         codetools_0.2-19     lubridate_1.9.3     
    [58] stringi_1.8.4        gtable_0.3.5         munsell_0.5.1       
    [61] tibble_3.2.1         pillar_1.9.0         htmltools_0.5.8.1   
    [64] ipred_0.9-15         IRkernel_1.3.2.9000  lava_1.8.0          
    [67] R6_2.5.1             evaluate_0.23        class_7.3-22        
    [70] Rcpp_1.0.12          uuid_1.2-0           nlme_3.1-162        
    [73] prodlim_2024.06.25   ModelMetrics_1.2.2.2 pkgconfig_2.0.3     

