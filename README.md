# CacaoGenomics
This repository contains scripts associated with the _Theobroma cacao_ genome sequencing project and related analyses.

## Usage

```
Usage
python3 exp_plots_tissue.py --genes <FILE> --exp <FILE> --out <FILE>

Mandatory:
--genes    STR   Candidate gene file.
--exp      STR   Count table.
--out      STR   Output folder
--samples  STR   Sample info file.

Optional:
--cutfac     INT Number of IQRs that determine outliers
--logscale   -   activates log scale [off]
--filteroff  -   Switches off outlier filter.

```

`--genes` specifies a text file containing the genes of interest. Each line lists one gene ID. These IDs need to match the IDs in the first column of the count table. The first column can be followed by additional columns with trivial names in a second column. Columns must be separated by tabs.

`--exp` specifies a count table with expression data with genes in rows and samples in columns.

`--out` specifies an output folder. If this folder does not exist, it will be created.

`--cutfac` specifies a multiple of the IQR that is used to classify samples as outliers.

`--logscale` activates the use of logarithmic transformation of all gene expression values.

`--filteroff` deactivates the outliers exclusion and shows all data points as a result.


## Reference

This repository.
