[![Github Release](https://img.shields.io/badge/release-v2.1-brightgreen)](https://github.com/3UTR/DaPars2)
[![python Release](https://img.shields.io/badge/python-3.8-brightgreen)](https://www.python.org/downloads/)
[![numpy Release](https://img.shields.io/badge/numpy-1.22-brightgreen)](https://numpy.org/)
[![scipy Release](https://img.shields.io/badge/numpy-1.80-brightgreen)](https://scipy.org/)
[![R Release](https://img.shields.io/badge/R-3.6.3-brightgreen)](https://cran.r-project.org/)


# DaPars2_LR (v1.0)

Dynamics analysis of Alternative PolyAdenylation from long-reads RNA-seq data, based on [DaPars2 (v2.1)](https://github.com/3UTR/DaPars2). 

## Introduction

The dynamic usage of the 3’untranslated region (3’UTR) resulting from alternative polyadenylation (APA) is emerging as a pervasive mechanism for regulating mRNA diversity, stability and translation. DaPars2_LR directly infers the dynamic alternative polyadenylation (APA) usage by comparing two samples. DaPars2_LR uses a different breakpoint detection method than DaPars2, described in our [publication](https://doi.org/10.1038/s41467-023-43387-9). It can infer the de novo proximal APA sites and returns long and short 3’UTR expression levels, associated with a FDR-adjusted Fisher test p-value. DaPars2_LR is the only tool for long-reads that uses coverage information (BedGraph), making it very efficient. 

![Flowchart](doc/Dapars_readme_COL6A1.png) 

## Workflow

![Flowchart](https://farm8.staticflickr.com/65535/51154541918_8a63879ed1_k.jpg)

## Usage

DaPars2_LR can be used just as DaPars2 (a [Full Documentations](https://github.com/3UTR/DaPars2/wiki) of the original DaPars2 can be found in Wiki page), just use `DaPars2_LR_Two_Samples_Multi_Chr.py` instead of `DaPars2_Multi_Samples_Multi_Chr.py`. **DaPars2_LR only accepts two-samples comparisons**.

DaPars2_LR uses a different breakpoint detection method, described in our [publication](https://doi.org/10.1038/s41467-023-43387-9).

Additionally, between step 1 and 2 of the [wiki](https://github.com/3UTR/DaPars2/wiki), DaPars2_LR offers to filter 3'UTR overlapping with a different gene, as they can create false positives:
```
python Dapars2_LR_Filter_Anno.py --input hg38_3UTR_annotation.bed --output hg38_3UTR_annotation_overlapping_filtered.bed
```

DaPars2_LR also offers to merge results from all Chr and compute FDR-adjusted p-values for each gene:
```
Dapars2_LR_merge.py -o Dapars2_LR_all_chr.txt
```

## Citation

*Please cite the following articles if you use DaPars2_LR in your research*
* Dondi, A. et al. Detection of isoforms and genomic alterations by high-throughput full-length single-cell RNA sequencing in ovarian cancer. NatCom (2023) [doi:110.1101/2022.12.12.520051.](https://doi.org/10.1038/s41467-023-43387-9)

*Please cite the following articles if you use DaPars2 in your research*:
* Feng X, Li L, Wagner EJ, Li W; TC3A: The Cancer 3′ UTR Atlas, Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D1027–D1030
* Li L^, Huang K^, Gao YP, Cui Y, Wang G, Nathan D, Li YM, Chen YE, Ji P, Peng F, William K, Wagner EJ, Li W. (2021) An atlas of alternative polyadenylation quantitative trait loci contributing to complex trait and disease heritability. Nature Genetics. doi: 10.1038/s41588-021-00864-5. 

## Contact

If you have any comments, suggestions, questions, bug reports, etc, feel free to contact Arthur Dondi (arthur.dondi@gmail.com). Please attach your command line and log messages if possible.
