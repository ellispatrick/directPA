### Direction analysis for pathways and kinases
========

#### Description
`directPA' is a package designed for analyzing (1) gene/protein expression data using microarray, RNA-seq, or proteomics data, and (2) kinase enrichment using phosphoproteomics data. It intergrates multiple experimental treatments and test the effect of each combination of treatments using a directional approach. 

#### Download and install
The relsease version can be downloaded from CRAN [link](http://cran.r-project.org/web/packages/directPA/);

Install the release version from CRAN with `install.packages("directPA")`

#### Examples
I. for kinase anlaysis on a phosphoproteomics dataset (Humphrey et al. Cell Metab., 2013)
```r
# load the phosphoproteomics dataset
data(HEK)

# load the kinase-substrate annoations
data(PhosphoSite)

# direction pathway analysis in 2-dimensional space. Implemented as rotating by degree
# (1) test combined effect of Torin1 and Rapamycin vs insul both on "down-regulation"
# (180 degree to original direction)
kst1 <- directPA(Tc=HEK, direction=pi, annotation=PhosphoSite.mouse)
kst1[order(unlist(kst1[,1])),][1:20,]

# (2) test combined effect of Torin1 and Rapamycin vs insul on "no change and down-regulation"
# (135 degree to the original direction)
kst2 <- directPA(Tc=HEK, direction=pi*3/4, annotation=PhosphoSite.mouse)
kst2[order(unlist(kst2[,1])),][1:20,]

# (3) test combined effect of Torin1 and Rapamycin vs insul on "down-regulation and no change"
# (225 degree to the original direction)
kst3 <- directPA(Tc=HEK, direction=pi*5/4, annotation=PhosphoSite.mouse)
kst3[order(unlist(kst3[,1])),][1:20,]
```
