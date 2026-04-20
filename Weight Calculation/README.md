## Workflow

1. Read annual raster layers for all indicators
2. Harmonize raster geometry to the same reference grid
3. Sample valid pixels for each year
4. Estimate within-group annual weights using CRITIC + entropy + MER fusion
5. Estimate between-group annual weights using CRITIC + entropy + MER fusion
6. Recover annual indicator weights
7. Estimate hierarchical fixed weights
8. Run spatial block bootstrap for uncertainty analysis


## Required packages

- terra
- dplyr
- purrr
- tibble
- readr
- tidyr
- yaml

## Usage

### 1) Copy and edit the config file

```r
file.copy("config.example.yml", "config.yml")
```

Then modify:
- input raster directories
- filename patterns
- output directory
- years
- grouping scheme
- bootstrap settings

### 2) Run the workflow

```bash
Rscript scripts/compute_hdci_weights.R config.yml
```

or in R:

```r
source("scripts/compute_hdci_weights.R")
run_hdci_weighting("config.yml")
```

## Core outputs

The script writes the following files to `output_dir`:

- `available_years.csv`
- `yearly_within_group_weights.csv`
- `yearly_between_group_weights.csv`
- `yearly_indicator_weights.csv`
- `fixed_between_group_weights.csv`
- `fixed_within_group_weights.csv`
- `final_fixed_weights.csv`
- `bootstrap_fixed_weights_summary.csv` (if enabled)

## Notes

1. All indicator rasters should already be transformed to the target impact scale (e.g., 0–100).
2. Annual raster files must share consistent naming rules defined in the config file.
3. The public release focuses on **weight estimation**, not final HDCI raster generation.
4. If needed, final HDCI rasters can be generated later by linearly combining annual indicator layers with `FinalFixedWeight_HierMER`.

## Citation

If you use this code in a publication, please cite the associated paper and describe the weighting workflow as:
“CRITIC + entropy + positive-constrained MER fusion with hierarchical fixed weighting.”
