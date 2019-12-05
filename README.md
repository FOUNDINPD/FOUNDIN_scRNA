# FOUNDIN_scRNA
Scripts/Code for analyzing single cell RNA-seq data from FOUNDIN PD.

There is a perl script to generate count matrix from Cell Ranger pipeline and a R script to integrate data from all the samples/cell lines.

SessionInfo for R

R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WebGestaltR_0.4.1      singleCellNet_0.1.0    reshape2_1.4.3         pheatmap_1.0.12        devtools_2.2.1         usethis_1.5.1         
 [7] gtools_3.8.1           stringi_1.4.3          corrr_0.3.2            doRNG_1.7.1            rngtools_1.4           pkgmaker_0.27         
[13] registry_0.5-1         foreach_1.4.4          GENIE3_1.6.0           SCENIC_1.1.2-2         readxl_1.3.1           jsonlite_1.6          
[19] Matrix_1.2-17          reshape_0.8.8          purrr_0.3.2            tximport_1.12.3        dplyr_0.8.3            reticulate_1.13.0-9000
[25] cowplot_0.9.4          ggplot2_3.2.0          Seurat_3.1.1          

loaded via a namespace (and not attached):
  [1] backports_1.1.4             plyr_1.8.4                  igraph_1.2.4.1              lazyeval_0.2.2              GSEABase_1.46.0            
  [6] splines_3.6.0               BiocParallel_1.18.0         listenv_0.7.0               GenomeInfoDb_1.20.0         digest_0.6.20              
 [11] htmltools_0.3.6             gdata_2.18.0                magrittr_1.5                memoise_1.1.0               doParallel_1.0.14          
 [16] cluster_2.0.8               ROCR_1.0-7                  remotes_2.1.0               readr_1.3.1                 globals_0.12.4             
 [21] annotate_1.62.0             RcppParallel_4.4.4          matrixStats_0.54.0          R.utils_2.9.0               prettyunits_1.0.2          
 [26] colorspace_1.4-1            blob_1.1.1                  ggrepel_0.8.1               callr_3.3.2                 crayon_1.3.4               
 [31] RCurl_1.95-4.12             graph_1.62.0                survival_2.44-1.1           zoo_1.8-6                   iterators_1.0.10           
 [36] ape_5.3                     glue_1.3.1                  gtable_0.3.0                zlibbioc_1.30.0             XVector_0.24.0             
 [41] leiden_0.3.1                DelayedArray_0.10.0         pkgbuild_1.0.5              future.apply_1.3.0          SingleCellExperiment_1.6.0 
 [46] apcluster_1.4.7             BiocGenerics_0.30.0         scales_1.0.0                DBI_1.0.0                   bibtex_0.4.2               
 [51] Rcpp_1.0.1                  metap_1.1                   viridisLite_0.3.0           xtable_1.8-4                bit_1.1-14                 
 [56] rsvd_1.0.1                  SDMTools_1.1-221.1          stats4_3.6.0                tsne_0.1-3                  htmlwidgets_1.3            
 [61] httr_1.4.0                  gplots_3.0.1.1              RColorBrewer_1.1-2          ellipsis_0.3.0              ica_1.0-2                  
 [66] pkgconfig_2.0.2             XML_3.98-1.20               R.methodsS3_1.7.1           uwot_0.1.4                  tidyselect_0.2.5           
 [71] rlang_0.4.0                 later_0.8.0                 AnnotationDbi_1.46.0        munsell_0.5.0               cellranger_1.1.0           
 [76] tools_3.6.0                 cli_1.1.0                   RSQLite_2.1.1               ggridges_0.5.1              stringr_1.4.0              
 [81] npsurv_0.4-0                processx_3.4.1              fs_1.3.1                    bit64_0.9-7                 fitdistrplus_1.0-14        
 [86] caTools_1.17.1.2            randomForest_4.6-14         RANN_2.6.1                  pbapply_1.4-0               future_1.14.0              
 [91] nlme_3.1-139                whisker_0.4                 mime_0.7                    R.oo_1.22.0                 compiler_3.6.0             
 [96] rstudioapi_0.10             plotly_4.9.0                png_0.1-7                   testthat_2.2.1              lsei_1.2-0                 
[101] tibble_2.1.3                ps_1.3.0                    desc_1.2.0                  lattice_0.20-38             pillar_1.4.2               
[106] Rdpack_0.11-0               lmtest_0.9-37               RcppAnnoy_0.0.12            data.table_1.12.2           bitops_1.0-6               
[111] irlba_2.3.3                 gbRd_0.4-11                 httpuv_1.5.1                AUCell_1.6.1                GenomicRanges_1.36.0       
[116] R6_2.4.0                    promises_1.0.1              KernSmooth_2.23-15          gridExtra_2.3               IRanges_2.18.1             
[121] sessioninfo_1.1.1           codetools_0.2-16            pkgload_1.0.2               MASS_7.3-51.4               assertthat_0.2.1           
[126] SummarizedExperiment_1.14.0 rprojroot_1.3-2             withr_2.1.2                 sctransform_0.2.0           S4Vectors_0.22.0           
[131] GenomeInfoDbData_1.2.1      hms_0.4.2                   parallel_3.6.0              grid_3.6.0                  tidyr_0.8.3                
[136] Rtsne_0.15                  Biobase_2.44.0              shiny_1.3.2                

