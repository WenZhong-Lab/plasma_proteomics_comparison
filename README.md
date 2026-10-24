# Comparative analysis of mass spectrometry- and affinity-based plasma proteomics platforms

This repository contains the R analysis code for a comparison of five plasma proteomics workflows: MS-OT480, MS-Astral, Olink-Explore, Olink-HT, and SomaScan. The study uses 30 plasma samples from 10 healthy donors, collected on days 0, 3, and 5.

The analyses examine protein detection and overlap, consistency across repeat samples, within- and between-donor variation, agreement between platforms, and protein signals associated with donor identity and sex. Scripts are grouped by the corresponding manuscript figures.

## Repository contents

| Path | Contents |
| --- | --- |
| `code/Figure_2.R`–`code/Figure_6.R` | Figure-specific analyses, summary tables, and plots |
| `code/setup.R` | Common paths and package loading |
| `code/functions.R`, `code/pdf_text_utils.R` | Shared analysis and figure-export helpers |
| `data/processed/` | Local processed assay inputs |
| `data/reference/` | Local reference annotations used by the analyses |
| `result/`, `plot/` | Generated analysis objects and figures |

## Working with the code

Run the figure scripts from the repository root in an R environment with the packages listed in `code/setup.R` and the additional packages loaded by individual figure scripts. `Figure_2.R` prepares detection summaries read by later scripts, so run it first when generating the figure outputs.

The supplied `.gitignore` keeps assay data, processed inputs, reference files, and generated results out of version control. The scripts use local input paths under `data/` and write outputs to `result/` and `plot/`.

## Study information

**Manuscript:** *Comparative analysis of mass spectrometry- and affinity-based plasma proteomics platforms*  
**Keywords:** plasma proteomics; mass spectrometry; Olink; SomaScan; cross-platform comparison
