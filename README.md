---
output:
  pdf_document: default
  html_document: default
---
# XVCF
<a href="https://github.com/rashidma/XVCF">XVCF</a> © 2025 by Ghaida Almuneef and Mamoon Rashid is licensed under <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nd.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">

XVCF (eXplore VCF) is a tool for visualization and quality control of genetic variation data.

## Motive
The premise behind developing "XVCF" shiny package is to pat the gap between
heap of genomic data available and their corresponding
visualization software/tool. 

## Solution
We addressed this need timely and provided easy to use,light-weight, GUI-based tool for the scientific community less versed in
programming or computer skills.

## Methods
We used "vcfR" and "maftools" R packages as back bones of "XVCF".

"XVCF" has two blocks operational. In the first block XVCF accepts a vcf file and generate some quality statistics like Read depth, Genotype, Genotype quality, Allele metrics etc. In the second module, it accepts ANNOVAR annoated vcf data in tabular format. The second module is basically for Cancer Genomic data and handled by underlying "maftools" R package. Below is the schematic diagram of XVCF depicting both the modules.

<img src="Fig1.png"> 

