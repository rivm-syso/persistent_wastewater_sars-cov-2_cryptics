# Persistent detection of two genetically and geographically distinct cryptic SARS-CoV-2 lineages in the Netherlands

Code for a publication currently under review.

## Description:

To be added.

## Installation:

The `environment` directory contains `def` files for Singularity/Apptainer containers which are used for the analysis (``` singularity build --fakeroot``environment``/name>.def environment/<name>.sif ``` ). The Snakemake pipeline should be run from a Conda environment generated using `environment/main.yaml`.

## Usage:

1.  Prepare the analysis with the `src/R/Prepare.R` script using the `R` Singularity/Apptainer container.

2.  Run the Snakemake pipeline using Singularity/Apptainer containers built from the environment files in `environment` .

3.  Perform the analysis with the `src/R/Analysis.R` script using the `R` Singularity/Apptainer container.

## Support:

If you encounter any problems, have questions, or would like to suggest improvements, please feel free to open an issue on this repository.

## Authorship:

### Code:

- [Auke Haver](https://orcid.org/0000-0002-6711-2205)

### Manuscript

- [Auke Haver](https://orcid.org/0000-0002-6711-2205)
- [Steff van Blokland](https://orcid.org/0009-0006-8100-0468)
- [Jaap T. van Dissel](https://orcid.org/0000-0002-3857-331X)
- [Jeroen F.J. Laros](https://orcid.org/0000-0002-8715-7371)
- [Willemijn J. Lodder](https://orcid.org/0009-0006-4795-5404)

## Software:

- tidyverse: <https://www.tidyverse.org/>
- ape: <https://cran.r-project.org/package=ape>
- ggtree: <https://bioconductor.org/packages/ggtree/>
- ggpubr: <https://cran.r-project.org/package=ggpubr>
- cowplot: <https://cran.r-project.org/package=cowplot>
- viridis: <https://cran.r-project.org/package=viridis>
- jsonlite: <https://cran.r-project.org/package=jsonlite>
- writexl: <https://cran.r-project.org/package=writexl>
- yaml: <https://cran.r-project.org/package=yaml>
- feather: <https://cran.r-project.org/package=feather>
- Snakemake: <https://snakemake.readthedocs.io/>
- Augur: <https://docs.nextstrain.org/projects/augur/>
- Freyja: <https://github.com/andersen-lab/Freyja>
- BiocManager: <https://cran.r-project.org/package=BiocManager>
