### Direction analysis for pathways and kinases
========

#### Description
`directPA' is a package designed to identify combinatorial effects of multiple treatments and/or perturbations on pathways and kinases profiled by microarray, RNA-seq, proteomics, or phosphoproteomics data.

#### Download and install
The relsease version can be downloaded from CRAN [link](https://cran.r-project.org/package=directPA);

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
kPA <- kinasePA(Tc=HEK, direction=pi, annotation=PhosphoSite.mouse)
kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(kPA$substrate.pvalues)[1:20]

# (2) test combined effect of Torin1 and Rapamycin vs insul on "no change and down-regulation"
# (135 degree to the original direction) 
kPA <- kinasePA(Tc=HEK, direction=pi*3/4, annotation=PhosphoSite.mouse)
kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]

# (3) test combined effect of Torin1 and Rapamycin vs insul on "down-regulation and no change"
# (225 degree to the original direction) 
kPA <- kinasePA(Tc=HEK, direction=pi*5/4, annotation=PhosphoSite.mouse)
kPA$kinase[order(unlist(kPA$kinase[,1])),][1:20,]
```
