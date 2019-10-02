# A comparison of short-read, long-read, and hybrid *de novo* genome assembly tools for the reconstruction of *Mycoplasma bovis* strains

We compared different assembly algorithms combined in several pipelines to reconstruct the genomes of *Mycoplasma bovis* strains.
Basically, we compared three assembly approaches: 
* only using short reads
* only using long reads
* using both read types

To sum it up, the use of short reads only lead to fragmented assemblies, but the sequences are rather accurate in relation to completely discovered genes (Snakefile and other files in folder illumina).
Long-read assembler can provide more contiguous results, but the higher error rate of MinION ONT data is a problem to get reliable results (Snakefile and other files in folder nanopore).
Thus, we noticed more indels, frame shifts or pseudogenes.
Contrary to our expectations, currently available de novo genome assembly tools with a hybrid option (e.g. SPAdes) did not perform well. <!--The fundamental principle of these tools, by building the assembly graph with short reads and improving it with long reads, cannot solve large repetitive regions.
While combining several methods with a focus on hybrid assembly (Unicycler) can improve the contiguity of assemblies, but long-read-only approaches are still performing better than short-read and hybrid approaches regarding this metric.-->
The second hybrid assembly approach combines the contiguity of the long-read assembler and ensures more accurate sequences using short reads in postprocessing (Snakefile and other files of the first and second approach in folder hybrid).
Finally, our assembly pipeline consists of four steps: short-read preprocessing, long-read preprocessing (demultiplexing and optionally filtering), long-read-only assembly, postprocessing with long and short reads.
Aftwerwards, QUAST is used to get some statistics to evaluate the assembly quality.
The correspondig Snakefile, config.yaml, create.py and conda environment files you can use to process you data you can find here.

How to use our pipeline:

0. install conda and snakemake
1. insert all paths and strains in the config.yaml
2. change something in the Snakefile
  * line 17: insert the name of one strain (to start the rule for demultiplexing right)
  * line 215: insert the path to the fastq file of the long reads
  * line 239 - rule rename_qcat: you have to insert the name of the strain corresponding to the used barcode, e.g. for strain 1234 barcode 1 was used, you have to insert 1234 in variable 'strain1' to replace STRAIN1 (if you want to use porechop, you have to uncomment these rules and update the rename rule for porechop)
  * line 335 and 403: change the path to the docker folder (docker is used for the Racon runs and docker has to have permissions on this folder, otherwise it would not work)
  * line 437: you have to insert all resulting assembly files as input to the quast rule to incorporate all of them in the QUAST statistics. Be aware that the number of input files is correct in the quast command in line 465.
3. run 'snakemake create' to create all folders needed in the pipeline steps
4. copy or link the files of the paired-end short reads to the folder /hybrid/raw_data/. The files should be named like this: STRAIN1_1.fastq.qz, STRAIN1_2.fastq.gz, ...
5. start the pipeline with 'snakemake -pr --cores &lt;threads&gt; --use-conda' (you can do a dry run of snakemake first to see, if all rules are going to be used with -prn instead of -pr)

