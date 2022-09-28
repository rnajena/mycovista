# TOOLNAME - a pipeline to assemble (highly repetitive) bacterial genomes

## Installation
In order to run TOOLNAME, you'll need to install:
* [conda](https://docs.conda.io/en/latest/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)

## Usage
To use TOOLNAME, edit the `config.yaml` file
* specify the pipeline you want to use: long or hybrid mode
* insert additional information:
  * path of the input and output folders
  * names of the bacterial strains

A `python` script helps you to generate several folders for the output. Just run:

``python scripts/create.py``

Afterwards, all required data should be linked in the `raw_data` folder. Please check this!

Now you can get started and assemble your bacterias. You can start TOOLNAME by using:

``snakemake -s <mode> -c <#threads> --use-conda``
* `mode - assembly mode file (long or hybrid)`

*Tip: You can use `-n` for a dry run to check if snakemake will start all required rules.*


## Reference list
Please cite the following tools when using our pipleine.

Short-read preprocessing:
* [fastp](https://github.com/OpenGene/fastp)
* [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)

Long-read preprocessing:
* [qcat](https://github.com/nanoporetech/qcat)
* [Filtlong](https://github.com/rrwick/Filtlong)

Long-read-only assembly:
* [Flye](https://github.com/fenderglass/Flye)

Assembly postprocessing:
* [Racon](https://github.com/isovic/racon)
* [minimap2](https://github.com/lh3/minimap2)
* [Medaka](https://github.com/nanoporetech/medaka)

Read quality evaluation:
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [NanoPlot](https://github.com/wdecoster/NanoPlot)

Assembly quality evaulation:
* [QUAST](https://github.com/ablab/quast)

Annotation:
* [Prokka](https://github.com/tseemann/prokka)
