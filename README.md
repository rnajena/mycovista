# A comparison of short-read, long-read, and hybrid *de novo* genome assembly tools for the reconstruction of *Mycoplasma bovis* strains

We compared different assembly algorithms combined in several pipelines to reconstruct the genomes of *Mycoplasma bovis* strains.
Basically, we compared three assembly approaches: 
* only using short reads
* only using long reads
* using both read types

To sum it up, the use of short reads only lead to fragmented assemblies, but the sequences are rather accurate in relation to completely discovered genes (short-read-only files in  'illumina' folder).
Long-read assembler can provide more contiguous results, but the higher error rate of MinION ONT data is a problem to get reliable results (long-read-only files in  'nanopore' folder).
Thus, we noticed more indels, frame shifts or pseudogenes.
Contrary to our expectations, currently available *de novo* genome assembly tools with a hybrid option (e.g. SPAdes) did not perform well.
<!--The fundamental principle of these tools, by building the assembly graph with short reads and improving it with long reads, cannot solve large repetitive regions.
While combining several methods with a focus on hybrid assembly (Unicycler) can improve the contiguity of assemblies, but long-read-only approaches are still performing better than short-read and hybrid approaches regarding this metric.-->
The second hybrid assembly approach combines the contiguity of the long-read assembler and ensures more accurate sequences using short reads in postprocessing (files of the first and second approach in 'hybrid' folder).
Finally, our assembly pipeline consists of four steps: short-read preprocessing, long-read preprocessing (demultiplexing and filtering), long-read-only assembly, postprocessing with long and short reads.
Aftwerwards, QUAST is used to get some statistics to evaluate the assembly quality and prokka generates annotation files using gene code 4.
The correspondig files (Snakefile, config.yaml, create.py conda environment files and scripts) you can find here.

How to use our pipeline:

0. Install conda and snakemake.
1. Download the reference genome (.fasta) and annotation (.gff) file.
2. Insert paths and strains in the config.yaml.
    * additional information is provided in the config.yaml file
3. Run `snakemake create` to create all folders needed during the pipeline steps.
4. Copy (`cp $SOURCE$ $DESTINATION$`) or link (`ln -s $SOURCE$ $DESTINATION$`) the raw data files to the (previous created) folder /raw_data/
    * paired-end short reads
      * The files should be named like this: STRAIN1_1.fastq.qz, STRAIN1_2.fastq.gz,...
    * long reads
      * The file should be named nanopore.fastq
5. Start the pipeline with `snakemake --use-conda --cores <threads>`.
    * You have to start snakemake while being in the git folder
    * You can start a dry run of snakemake first to see, if all rules are going to be used with `-prn`
    * If you want to see additional information (e.g. why this rule is started,..), you can add `-pr` (print reason) to the snakemake command
