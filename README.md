# analysePeaks
A package to compare your peaks for ChIP/CUT&amp;RUN/ATAC... between them
This package will compare peaks (usually narrow peaks). It can be comparing replicates or comparing different techniques. Compare with peaks from another experiment. Compare to TSS, to binding motifs...

## Installation
This package depends on grDevices, graphics, utils, usefulLDfunctions, usefulLDfunctionsGR, GenomicRanges, BiocGenerics, stats, pheatmap, reshape.
If usefulLDfunctions, and/or usefulLDfunctionsGR are not installed on your R see Installation of dependencies section.

The easiest way to install analysePeaks is using devtools::install_github() from R:
```
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/analysePeaks")
# Or, if you want to have the tutorial (this needs knitr, see Installation of dependencies section):
devtools::install_github("lldelisle/analysePeaks", build_vignettes = TRUE)
```

## Tutorial
A vignette for this package have been written.
To see it, you need to run:
```
vignette("getting-started-with-analysePeaks")
```
If you get the message: `vignette ‘getting-started-with-analysePeaks’ not found`. You need to reinstall the package with the vignette:
```
devtools::install_github("lldelisle/analysePeaks", build_vignettes = TRUE, force = TRUE)
```



## Installation of dependencies
usefulLDfunctions is another R package I wrote which is on github too. So to install it:
```
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
```
As the installation of Bioconductor package depends on the R version you have, I recommend you to use the function `safelyLoadAPackageInCRANorBioconductor` from the usefulLDfunctions package:
```
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("BiocGenerics")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("combinat")
safelyLoadAPackageInCRANorBioconductor("knitr")
```
usefulLDfunctionsGR is another R package I wrote which is on github too. So to install it:
```
devtools::install_github("lldelisle/usefulLDfunctionsGR")
```

## Issues
If you have issues, use the "Issues" section in GitHub or send an email to lucille.delisle\@epfl.ch
