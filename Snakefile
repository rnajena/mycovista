# de novo hybrid assembly pipeline

# files and directories
configfile: "config.yaml"

PU = '1P 1U 2P 2U'.split() # paired_unpaired
END = '1 2'.split() # paired_end

# barcode_strain = dict(config["strains"])

# strains = list(barcode_strain.keys())
strains = list(config["strains"])
# barcodes = list(barcode_strain.values())
# for elem in strains[:]:
# 	if elem[0:6] == 'strain':
# 		help = strains.index(elem)
# 		strains.remove(elem)
# 		barcodes.pop(help)
# 
# refgenome = config["refgenome"][0]
# refannotation = config["refannotation"][0]
refdatabase = config["refdatabase"][0]


rule all:
	input:
		# qcat
		# expand("{path}/preprocessing/qcat/none.fastq", path = config["path"]),
		# expand("{path}/preprocessing/qcat/{strain}_qcat.fastq", path = config["path"], strain = strains),
		# 
		# filtlong
		expand("{path}/preprocessing/{demultiplex}/{strain}_{demultiplex}_filtered.fastq.gz", path = config["path"],  strain = strains, demultiplex = config["demultiplexing"]),
		# 
		# nanoplot
		expand("{path}/quality/nanoplot/{strain}/{strain}_{demultiplex}_NanoPlot-report.html", path = config["path"],  strain = strains, demultiplex = config["demultiplexing"]),
		# 
		# preprocess short reads
		expand("{path}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.html", path = config["path"], strain = strains, paired_end = END),
		expand("{path}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.html", path = config["path"], strain = strains, paired_unpaired = PU),
		expand("{path}/preprocessing/illumina/{strain}_unique.fastq", path = config["path"], strain = strains),
		# 
		# assembly
		expand("{path}/assembly/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}.fasta", path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# 
		# polishing - 4x Racon long -> medaka -> 4x Racon short
		expand("{path}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_long4.fasta", path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		expand("{path}/postprocessing/{strain}_{demultiplex}_{assembler}/consensus.fasta", path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		expand("{path}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_final.fasta", path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		
		# quast
		expand("{path}/quality/quast/report.html", path = config["path"], strain = strains),
		# 
		# prokka
		# expand("{path}/annotation/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}.gff", path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		#
		# ideel
		# expand('{path}/ideel/hists/{strain}_{demultiplex}_{assembler}_final.pdf',  path = config["path"], strain = strains, demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# 
		# 
		# quast
		# expand("{path}/quast/report.html", path = config["path"], strain = strains),
		# 
		# prokka
		expand("{path}/annotation/{strain}.gff", path = config["path"], strain = strains),
		#
		# ideel
		# expand('{path}/ideel/hists/{strain}.pdf',  path = config["path"], strain = strains)


# create folders for following steps
rule create:
	shell:
		'python scripts/create.py'

# quality check of the short raw reads - FastQC
rule fastqc:
	input:
		'{path}/raw_data/{strain}_{paired_end}.fastq.gz'
	output:
		html = '{path}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.html',
		zip = '{path}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.zip'
	conda:
		'envs/read_quality.yaml'
	params:
		outputdir = '{path}/quality/fastqc/{strain}/'
	threads: 8
	shell:
		'fastqc -t {threads} {input} -o {params.outputdir}'


# first preprocessing of the short reads - fastp
rule fastp:
	input:
		fw = '{path}/raw_data/{strain}_1.fastq.gz',
		rv = '{path}/raw_data/{strain}_2.fastq.gz'
	output:
		f = '{path}/preprocessing/illumina/fastp/{strain}_1_fastp.fastq.gz',
		r = '{path}/preprocessing/illumina/fastp/{strain}_2_fastp.fastq.gz'
	conda:
		'envs/preprocessing.yaml'
	params:
		outputdir = '{path}/preprocessing/illumina/fastp/',
		html = '{strain}_fastp.html',
		json = '{strain}_fastp.json'
	threads: 8
	shell:
		'fastp -w {threads} -i {input.fw} -I {input.rv} -o {output.f} -O {output.r}' # --json {params.outputdir}{params.json} --html {params.outputdir}{params.html}'


# second preprocessing of the short reads - Trimmomatic
rule trimmomatic:
	input:
		fw = rules.fastp.output.f,
		rv = rules.fastp.output.r
	output:
		fP = '{path}/preprocessing/illumina/trimmomatic/{strain}_1P.fastq.gz',
		fU = '{path}/preprocessing/illumina/trimmomatic/{strain}_1U.fastq.gz',
		rP = '{path}/preprocessing/illumina/trimmomatic/{strain}_2P.fastq.gz',
		rU = '{path}/preprocessing/illumina/trimmomatic/{strain}_2U.fastq.gz'
	conda:
		'envs/preprocessing.yaml'
	threads: 8
	shell:
		'trimmomatic PE -phred33 -threads {threads} {input.fw} {input.rv} {output.fP} {output.fU} {output.rP} {output.rU} SLIDINGWINDOW:4:28 MINLEN:20'
			# PE					for paired end reads
			# -phred33				Illumina 1.9 uses phred +33
			# SLIDINGWINDOW:4:28	quality score = 28, because the trimming can be done more restrictive due to the high coverage    window size = 4 to prevent to stop the trimming at a local score maximum within one read
			# MINLEN:20				to get less random matches during mapping


# concatenate short reads for further postprocessing of the assembly
rule concat_short:
	input:
		fwP = rules.trimmomatic.output.fP,
		fwU = rules.trimmomatic.output.fU,
		rvP = rules.trimmomatic.output.rP,
		rvU = rules.trimmomatic.output.rU
	output:
		all = '{path}/preprocessing/illumina/{strain}_all_short.fastq.gz'
	shell:
		'cat {input.fwP} {input.fwU} {input.rvP} {input.rvU} > {output.all}'

# unzip short read files for further processing
rule gunzip_short:
	input:
		all = rules.concat_short.output.all
	output:
		gunzip = '{path}/preprocessing/illumina/{strain}_all_short.fastq'
	shell:
		'gunzip {input.all}'


# concatenate read ID to get unique read IDs for further postprocessing of the assembly
rule unique_readID:
	input:
		all = rules.gunzip_short.output.gunzip
	output:
		unique = '{path}/preprocessing/illumina/{strain}_unique.fastq'
	shell:
		"sed 's/ 2:N:0/:2:N:0/g' {input.all} | sed 's/ 1:N:0/:1:N:0/g' > {output.unique}"


# quality check of the preprocessed short reads - FastQC
rule fastqc_preprocessing:
	input:
		'{path}/preprocessing/illumina/trimmomatic/{strain}_{paired_unpaired}.fastq.gz'
	output:
		html='{path}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.html',
		zip='{path}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.zip'
	conda:
		'envs/read_quality.yaml'
	params:
		outputdir = '{path}/quality/fastqc/{strain}/'
	threads: 8
	shell:
		'fastqc {input} -t {threads} -o {params.outputdir}'


# demultiplexing of the long reads - qcat
rule qcat:
	input:
		'{path}/raw_data/nanopore.fastq'
	output:
		'{path}/preprocessing/qcat/none.fastq'
		# ...
	conda:
		'envs/qcat.yaml'
	params:
		outputdir = '{path}/preprocessing/qcat/'
	threads: 16 # default: 8
	shell:
		'qcat -t {threads} -k NBD104/NBD114 --trim -f {input} -b {params.outputdir}'


# rename qcat output
def get_input(strain):
    return '/preprocessing/qcat/' + barcodes[strains.index(strain)] + '.fastq'

rule rename_qcat:
    input:
        '{path}/preprocessing/qcat/none.fastq'
    output:
        '{path}/preprocessing/qcat/{strain}_qcat.fastq'
    params:
        strain = '{strain}',
        path = '{path}'
    run:
        import os
        command = 'mv ' + str(params.path) + get_input(str(params.strain)) + ' ' + str(output)
        print(command + '\n')
        os.system(command)


# filter demultiplexed long reads - Filtlong
rule filtlong:
	input:
		reads = '{path}/raw_data/{strain}_{demultiplex}.fastq'
		# reads = '{path}/preprocessing/{demultiplex}/{strain}_{demultiplex}.fastq'
	output:
		filtered = '{path}/preprocessing/{demultiplex}/{strain}_{demultiplex}_filtered.fastq.gz'
	conda:
		'envs/preprocessing.yaml'
	threads: 8
	shell:
		'filtlong --min_length 1000 {input.reads} | gzip > {output.filtered}'


# quality check of the preprocessed long reads - NanoPlot
rule nanoplot:
	input:
		rules.filtlong.output.filtered
	output:
		'{path}/quality/nanoplot/{strain}/{strain}_{demultiplex}_NanoPlot-report.html'
	conda:
		'envs/read_quality.yaml'
	params:
		outputdir = '{path}/quality/nanoplot/{strain}/',
		prefix = '{strain}_{demultiplex}_'
	threads: 8
	shell:
		'NanoPlot -t {threads} --minlength 1000 --fastq {input} -o {params.outputdir} -p {params.prefix} --title {params.prefix}minlength1000'


# long-read-only assembly - Flye
rule flye:
	input:
		rules.filtlong.output.filtered
	output:
		contigs = '{path}/assembly/{strain}_{demultiplex}_flye/assembly.fasta'
	conda:
		'envs/assembly.yaml'
	params:
		outputdir = '{path}/assembly/{strain}_{demultiplex}_flye/',
	threads: 16
	shell:
		'flye --asm-coverage 50 --iterations 2 --plasmids --nano-raw {input} -o {params.outputdir} -t {threads} -g 1000000'
        # --asm-coverage 50         coverage: 50X
        # -g 1000000                genome size 1 million
		# --iterations 2			two runs of polishing
		# --plasmids				for circular genomes and plasmids


# rename Flye output
rule rename_flye:
	input:
		contigs = rules.flye.output.contigs
	output:
		flye = '{path}/assembly/{strain}_{demultiplex}_flye/{strain}_{demultiplex}_flye.fasta'
	shell:
		'mv {input.contigs} {output.flye}'


# assembly postprocessing (polishing) using long reads - Racon			including minimap2 for the mapping inbetween
rule minimap2_racon_long:
	input:
		assembly = rules.rename_flye.output.flye,
		reads = rules.filtlong.output.filtered
	output:
		out = '{path}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_long4.fasta'
	conda:
		'envs/postprocessing.yaml'
	params:
		strain = '{strain}',
		demultiplex = '{demultiplex}',
		assembler = '{assembler}',
		path = '{path}'
	threads: 16
	script:
		'scripts/racon_long.py'


# assembly postprocessing (polishing) using long reads - medaka
rule medaka:
	input:
		reads = rules.filtlong.output.filtered,
		racon_out = rules.minimap2_racon_long.output.out
	output:
		medaka = '{path}/postprocessing/{strain}_{demultiplex}_{assembler}/consensus.fasta'
	conda:
		'envs/postprocessing.yaml'
	params:
		outputdir = '{path}/postprocessing/{strain}_{demultiplex}_{assembler}/'
	threads: 16
	shell:
		'medaka_consensus -i {input.reads} -d {input.racon_out} -o {params.outputdir} -t {threads} -m r941_min_high_g344'


# assembly postprocessing (polishing) using short reads - Racon			including minimap2 for the mapping inbetween
rule minimap2_racon_short:
	input:
		assembly = rules.medaka.output.medaka,
		reads = rules.unique_readID.output.unique
	output:
		out = '{path}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_short4.fasta'
	conda:
		'envs/postprocessing.yaml'
	params:
		strain = '{strain}',
		demultiplex = '{demultiplex}',
		assembler = '{assembler}',
		path = '{path}'
	threads: 16
	script:
		'scripts/racon_short.py'


# rename final assembly
rule final:
	input:
		rules.minimap2_racon_short.output.out
	output:
		'{path}/final_assemblies/{strain}.fasta'
	shell:
		'mv {input} {output}'


# assembly annotation - Prokka	
rule prokka:
	input:
		assembly = '{path}/final_assemblies/{strain}.fasta'
	output:
		gff = '{path}/annotation/{strain}.gff'
	conda:
		'envs/annotation.yaml'
	params:
		outputdir = '{path}/annotation/',
		prefix = '{strain}'
	threads: 16
	shell:
		'prokka --cpus {threads} --gcode 4 --force --outdir {params.outputdir} --prefix {params.prefix} {input.assembly}'


# assembly statistics (of all generated assemblies) - QUAST
def get_quast_in():
	out = ''
	for i in strains:
		out += config["path"][0] + '/final_assemblies/' + i + '.fasta '
	return out[0:len(out)-1]
quast_input = get_quast_in().split(" ")

rule quast:
	input:
		quast_input
	output:
		report = '{path}/quality/quast/report.html'
	conda:
		'envs/assembly_quality.yaml'
	params:
		outputdir = '{path}/quality/quast/',
		i = get_quast_in()
	threads: 16
	shell:
		'quast {params.i} -o {params.outputdir} -t {threads}' # -r ' + refgenome + ' -g ' + refannotation


# ideel plot
# Quick test for interrupted ORFs in bacterial/microbial genomes - https://github.com/phiweger/ideel
rule prodigal:
    input:
        '{path}/final_assemblies/{strain}.fasta'
    output:
        temp('{path}/ideel/proteins/{strain}.faa')
    conda:
        'envs/ideel.yaml'
    shell:
        'prodigal -g 4 -a {output} -q -i {input}'

rule diamond:
    input:
        '{path}/ideel/proteins/{strain}.faa'
    output:
        temp('{path}/ideel/lengths/{strain}.data')
    conda:
        'envs/ideel.yaml'
    threads: 8
    params:
        of = '6 qlen slen'
    shell:
        'diamond blastp --threads {threads} --max-target-seqs 1 --db ' + refdatabase + ' --query {input} --outfmt {params.of} --out {output}'

rule hist:
    input:
        '{path}/ideel/lengths/{strain}.data'
    output:
        '{path}/ideel/hists/{strain}.pdf'
    conda:
        'envs/ideel.yaml'
    script:
        'scripts/hist.R'
    # http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts
