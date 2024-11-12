# CacaoGenomics
scripts associated with the Theobroma cacao genome sequencing project

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


## Reference

This repository.
