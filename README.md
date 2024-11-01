# HGEN677 Presentation

### Dependence and Independence in Polygenic Scores
[![Licence](https://img.shields.io/github/license/Ileriayo/markdown-badges?style=for-the-badge)](./LICENSE)
#### General Information

- **Instructor**: Celia M Greenwood  
- **Presenter**: Arsham Mikaeili Namini  
- **Date**: October 29, 2024  
- **Department**: Human Genetics, Faculty of Medicine and Health Sciences, McGill University  

---

### Outline
#### Papers

1. **Pavlidis and Noble (2001)**  
   _Analysis of strain and regional variation in gene expression in mouse brain._  
   Published in **Genome Biology**, Volume 2: research0042.1  
   [Read Article](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2001-2-10-research0042)

2. **Aarts et al. (2014)**  
   _A solution to dependency: using multilevel analysis to accommodate nested data._  
   Published in **Nature Neuroscience**, Volume 17, Pages 491–496 <br/>
   [Read Article](https://www.nature.com/articles/nn.3648)

#### Focus

- **Figures 1 and 2**: In-depth explanation of results in Aarts et al. (2014).
- **Methodologies**: Overview of linear mixed model approaches, particularly GCTA and similar tools.

---

### Materials

This presentation uses a previously published set of gene expression microarray data from six brain regions in two mouse strains. In the [previous analysis](https://www.pnas.org/doi/abs/10.1073/pnas.97.20.11038), 24 genes showed strain-based expression differences, and around 240 genes exhibited regional expression differences.

The data, analysis results, and methodology are accessible from [this page](https://home.pavlab.msl.ubc.ca/mousebrain/).

To download the raw microarray data directly:

```bash
wget https://home.pavlab.msl.ubc.ca/mousebrain/Alldatanomef.txt
```
***
### Outcome
Using a mixed-effects linear model instead of conventional methods like t-tests and ANOVA, we reduced Type I error inflation and increased statistical power by assuming biological replication as a random effect.

#### Presentation slides are available on Google Slides:
[![Google Slides Badge](https://img.shields.io/badge/Google%20Slides-Presentation-yellow?style=for-the-badge&logo=google-slides&logoColor=white)](https://docs.google.com/presentation/d/1RqewBCk7_UbwAoM_LwZ-ef-t1a7AZ8HVd-8APOeXPRo/edit?usp=sharing)
***
### Reproducibility

To reproduce this analysis in R, follow the steps below:

1. Clone the Repository

Start by cloning the repository to your local machine:
```bash
git clone git@github.com:Arshammik/HGEN677_presentation.git
cd HGEN677-presentation
```

2. Install Required Packages

Open R or RStudio and set the working directory to the cloned repository folder. Use the following command in R to install required packages:

```r

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, ggplot2, dplyr, repr, pheatmap, irr, lmerTest, 
               gridExtra, pROC, countsplit, caret)


```
