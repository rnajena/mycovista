# A comparison of short-read, long-read, and hybrid *de novo* genome assembly tools for the reconstruction of *Mycoplasma bovis* strains

We compared different assembly algorithms combined in several pipelines to reconstruct the genomes of *Mycoplasma bovis* strains.
Basically, we compared three assembly approaches: 
* only using short reads
* only using long reads
* using both read types

To sum it up, the use of short reads only lead to fragmented assemblies, but the sequences are rather accurate in relation to completely discovered genes.
Long-read assembler can provide more contiguous results, but the higher error rate of MinION ONT data is a problem to get reliable results.
Thus, we noticed more indels, frame shifts or pseudogenes.
Contrary to our expectations, currently available *de novo* genome assembly tools with a hybrid option (e.g. SPAdes) did not perform well.
<!--The fundamental principle of these tools, by building the assembly graph with short reads and improving it with long reads, cannot solve large repetitive regions.
While combining several methods with a focus on hybrid assembly (Unicycler) can improve the contiguity of assemblies, but long-read-only approaches are still performing better than short-read and hybrid approaches regarding this metric.-->
The second hybrid assembly approach combines the contiguity of the long-read assembler and ensures more accurate sequences using short reads in postprocessing.
Finally, our assembly pipeline consists of four steps: short-read preprocessing (filtering and trimming), long-read preprocessing (demultiplexing and filtering), long-read-only assembly, postprocessing with long and short reads.
Aftwerwards, QUAST is used to generate some statistics to evaluate the assembly quality and prokka generates annotation files using gene code 4.
Moreover, this pipeline generates ideel plots of the final assemblies.


## How to use our pipeline:

0. Install conda and snakemake.
1. Download some reference files.
    * You need to download a reference genome (.fasta) and reference annotation (.gff) file.
    * Additional you will need to create a DIAMOND database as reference for the ideel plots. Therefore, you can use th following command: `diamond makedb --in database.fasta(.gz) --db database.dmnd`
2. Insert paths and strains in the config.yaml.
    * Additional information of the paths is provided in the config.yaml file.
3. Run `snakemake create` to create all folders needed during the pipeline steps.
4. Copy (`cp $SOURCE$ $DESTINATION$`) or link (`ln -s $SOURCE$ $DESTINATION$`) the raw data files to the (previous created) folder /raw_data/
    * paired-end short reads
      * The files should be named like: STRAIN1_1.fastq.qz, STRAIN1_2.fastq.gz,...
    * long reads
      * The file should be named: nanopore.fastq
5. Start the pipeline with `snakemake --use-conda --cores <threads>`.
    * You have to start snakemake while being in the git folder.
    * You can start a dry run of snakemake first to see, if all rules are going to be used with `-prn`.
    * If you want to see additional information (e.g. why this rule is started,..), you can add `-pr` (print reason) to the snakemake command.


## Reference list
Please cite the following tools when using our pipleine.

General tools:
* [conda](https://anaconda.org/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)

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

Visualization:
* [ideel](https://github.com/phiweger/ideel)
  * [Prodigal](https://github.com/hyattpd/Prodigal)
  * [Diamond](https://github.com/bbuchfink/diamond)
  * R (including libraries readr and ggplot2)
