# hybrid assembly pipeline

# files and directories
configfile: "config.yaml"

main_folder = 'hybrid'

PU = '1P 1U 2P 2U'.split() # paired_unpaired
END = '1 2'.split() # paired_end

rule all:
	input:
		# porechop
		# expand("{path}/{main_folder}/preprocessing/porechop/14DD0147_porechop.fastq", path = config["path"], main_folder = main_folder),
		# qcat
		# expand("{path}/{main_folder}/preprocessing/qcat/{strain}_qcat.fastq", path = config["path"], main_folder = main_folder, strain = config["strains"].values()),
        # 
		# NanoPlot
		# expand("{path}/{main_folder}/quality/nanoplot/{strain}/{strain}_{demultiplex}NanoPlot-report.html", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"]),
		# # 
		# # preprocess short reads
		# expand("{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.html", path = config["path"], main_folder = main_folder, strain = config["strains"], paired_end = END),
		# expand("{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.html", path = config["path"], main_folder = main_folder, strain = config["strains"], paired_unpaired = PU),
		# expand("{path}/{main_folder}/preprocessing/{strain}_unique.fastq", path = config["path"], main_folder = main_folder, strain = config["strains"]),
		# # 
		# # assembly
		# expand("{path}/{main_folder}/assembly/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}.fasta", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# # 
		# # polishing - 4x Racon long -> medaka -> 4x Racon short
		# expand("{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_long4.fasta", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# expand("{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/consensus.fasta", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# expand("{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_short4.fasta", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"], assembler = config["assembly"]),
		# # 
		# # QUAST
		# expand("{path}/{main_folder}/quality/quast/report.html", path = config["path"], main_folder = main_folder, strain = config["strains"]),
		# # 
		# # Prokka
		# expand("{path}/{main_folder}/prokka/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_short4.gff", path = config["path"], main_folder = main_folder, strain = config["strains"], demultiplex = config["demultiplexing"], assembler = config["assembly"])


# create folders for the following steps
rule create:
	shell:
		'python create.py'

# quality check of all short raw reads
rule fastqc:
	input:
		'{path}/{main_folder}/raw_data/{strain}_{paired_end}.fastq.gz'
	output:
		html = '{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.html',
		zip = '{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_end}_fastqc.zip'
	conda:
		'read_quality.yml'
	params:
		outputdir = '{path}/{main_folder}/quality/fastqc/{strain}/'
	threads: 8
	shell:
		'fastqc -t {threads} {input} -o {params.outputdir}'

# first preprocessing of the short reads - fastp
rule fastp:
	input:
		forward = '{path}/{main_folder}/raw_data/{strain}_1.fastq.gz',
		reverse = '{path}/{main_folder}/raw_data/{strain}_2.fastq.gz'
	output:
		forward = '{path}/{main_folder}/preprocessing/fastp/{strain}_1_fastp.fastq.gz',
		reverse = '{path}/{main_folder}/preprocessing/fastp/{strain}_2_fastp.fastq.gz'
	conda:
		'preprocessing.yml'
	params:
		outputdir = '{path}/{main_folder}/preprocessing/fastp/',
		html = '{strain}_fastp.html',
		json = '{strain}_fastp.json'
	threads: 8
	shell:
		'fastp -w {threads} -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} --json {params.outputdir}{params.json} --html {params.outputdir}{params.html}'

# second preprocessing of the short reads - trimmomatic
rule trimmomatic:
	input:
		forward = rules.fastp.output.forward,
		reverse = rules.fastp.output.reverse
	output:
		forwardPaired = '{path}/{main_folder}/preprocessing/trimmomatic/{strain}_1P.fastq.gz',
		forwardUnpaired = '{path}/{main_folder}/preprocessing/trimmomatic/{strain}_1U.fastq.gz',
		reversePaired = '{path}/{main_folder}/preprocessing/trimmomatic/{strain}_2P.fastq.gz',
		reverseUnpaired = '{path}/{main_folder}/preprocessing/trimmomatic/{strain}_2U.fastq.gz'
	conda:
		'preprocessing.yml'
	threads: 8
	shell:
		'trimmomatic PE -phred33 -threads {threads} {input.forward} {input.reverse} {output.forwardPaired} {output.forwardUnpaired} {output.reversePaired} {output.reverseUnpaired} SLIDINGWINDOW:4:28 MINLEN:20'
			# PE					for paired end reads
			# -phred33				Illumina 1.9 uses phred +33
			# SLIDINGWINDOW:4:28	quality score = 28, because the trimming can be done more restrictive due to the high coverage    window size = 4 to prevent to stop the trimming at a local score maximum within one read
			# MINLEN:20				to get less random matches during mapping

# concatenate all short reads separated by strain for postprocessing of the long reads assemblies later
rule concat_short:
	input:
		forwardPaired = rules.trimmomatic.output.forwardPaired,
		forwardUnpaired = rules.trimmomatic.output.forwardUnpaired,
		reversePaired = rules.trimmomatic.output.reversePaired,
		reverseUnpaired = rules.trimmomatic.output.reverseUnpaired
	output:
		all = '{path}/{main_folder}/preprocessing/{strain}_all_short.fastq.gz'
	shell:
		'cat {input.forwardPaired} {input.forwardUnpaired} {input.reversePaired} {input.reverseUnpaired} > {output.all}'

rule gunzip_short:
	input:
		all = rules.concat_short.output.all
	output:
		gunzip = '{path}/{main_folder}/preprocessing/{strain}_all_short.fastq'
	shell:
		'gunzip {input.all}'

# concatenate the read ID in one line to get unique read IDs for postprocessing of the long reads assemblies later
rule unique_readID:
	input:
		all = rules.gunzip_short.output.gunzip
	output:
		unique = '{path}/{main_folder}/preprocessing/{strain}_unique.fastq'
	shell:
		"sed 's/ 2:N:0/:2:N:0/g' {input.all} | sed 's/ 1:N:0/:1:N:0/g' > {output.unique}"

# quality check of all preprocessed short reads
rule fastqc_preprocessing:
	input:
		'{path}/{main_folder}/preprocessing/trimmomatic/{strain}_{paired_unpaired}.fastq.gz'
	output:
		html='{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.html',
		zip='{path}/{main_folder}/quality/fastqc/{strain}/{strain}_{paired_unpaired}_fastqc.zip'
	conda:
		'read_quality.yml'
	params:
		outputdir = '{path}/{main_folder}/quality/fastqc/{strain}/'
	threads: 8
	shell:
		'fastqc {input} -t {threads} -o {params.outputdir}'


long_path = config["long_path"].keys()
# long_path = long_path[2:len(long_path)-8]

# demultiplexing the long reads with qcat
rule qcat:
	input:
        '/data/fass2/reads/ont/celia/mycoplasma_doubice_2018/called/2019-06-20/myco_doubice_2019-06-20_2019-06-20_12-29-43_16193.fastq'
	output:
		BC01 = '{path}/{main_folder}/preprocessing/qcat/barcode01.fastq',
		BC02 = '{path}/{main_folder}/preprocessing/qcat/barcode02.fastq',
		BC03 = '{path}/{main_folder}/preprocessing/qcat/barcode03.fastq',
		BC04 = '{path}/{main_folder}/preprocessing/qcat/barcode04.fastq',
		BC05 = '{path}/{main_folder}/preprocessing/qcat/barcode05.fastq',
		BC06 = '{path}/{main_folder}/preprocessing/qcat/barcode06.fastq',
		BC07 = '{path}/{main_folder}/preprocessing/qcat/barcode07.fastq',
		BC08 = '{path}/{main_folder}/preprocessing/qcat/barcode08.fastq',
		BC09 = '{path}/{main_folder}/preprocessing/qcat/barcode09.fastq',
		BC10 = '{path}/{main_folder}/preprocessing/qcat/barcode10.fastq',
		BC11 = '{path}/{main_folder}/preprocessing/qcat/barcode11.fastq',
		BC12 = '{path}/{main_folder}/preprocessing/qcat/barcode12.fastq'
	conda:
		'qcat.yml'
	params:
		outputdir = '{path}/{main_folder}/preprocessing/qcat/'
	threads: 32 #default: 8
	shell:
		'qcat -t {threads} -k NBD103/NBD104 --trim -f {input} -b {params.outputdir}'

# rename qcat output
rule rename_qcat:
	input:
		BC01 = rules.qcat.output.BC01,
		BC02 = rules.qcat.output.BC02,
		BC03 = rules.qcat.output.BC03,
		BC04 = rules.qcat.output.BC04,
		BC05 = rules.qcat.output.BC05,
		BC06 = rules.qcat.output.BC06,
		BC07 = rules.qcat.output.BC07,
		BC08 = rules.qcat.output.BC08,
		BC09 = rules.qcat.output.BC09,
		BC10 = rules.qcat.output.BC10,
		BC11 = rules.qcat.output.BC11,
		BC12 = rules.qcat.output.BC12
	output:
		strain1 = '{path}/{main_folder}/preprocessing/qcat/STRAIN1_qcat.fastq',
		strain2 = '{path}/{main_folder}/preprocessing/qcat/STRAIN2_qcat.fastq',
		strain3 = '{path}/{main_folder}/preprocessing/qcat/STRAIN3_qcat.fastq',
		strain4 = '{path}/{main_folder}/preprocessing/qcat/STRAIN4_qcat.fastq',
		strain5 = '{path}/{main_folder}/preprocessing/qcat/STRAIN5_qcat.fastq',
		strain6 = '{path}/{main_folder}/preprocessing/qcat/STRAIN6_qcat.fastq',
		strain7 = '{path}/{main_folder}/preprocessing/qcat/STRAIN7_qcat.fastq',
		strain8 = '{path}/{main_folder}/preprocessing/qcat/STRAIN8_qcat.fastq',
		strain9 = '{path}/{main_folder}/preprocessing/qcat/STRAIN9_qcat.fastq',
		strain10 = '{path}/{main_folder}/preprocessing/qcat/STRAIN10_qcat.fastq',
		strain11 = '{path}/{main_folder}/preprocessing/qcat/STRAIN11_qcat.fastq',
		strain12 = '{path}/{main_folder}/preprocessing/qcat/STRAIN12_qcat.fastq'
	shell:
		'mv {input.BC01} {output.strain1} &&'
		'mv {input.BC02} {output.strain2} &&'
		'mv {input.BC03} {output.strain3} &&'
		'mv {input.BC04} {output.strain4} &&'
		'mv {input.BC05} {output.strain5} &&'
		'mv {input.BC06} {output.strain6} &&'
		'mv {input.BC07} {output.strain7} &&'
		'mv {input.BC08} {output.strain8} &&'
		'mv {input.BC09} {output.strain9} &&'
		'mv {input.BC10} {output.strain10} &&'
		'mv {input.BC11} {output.strain11} &&'
		'mv {input.BC12} {output.strain12}'


# filter demultiplexed long reads with Filtlong
rule filtlong:
	input:
		reads = '{path}/{main_folder}/preprocessing/{demultiplex}/{strain}_{demultiplex}.fastq'
	output:
		'{path}/{main_folder}/preprocessing/{demultiplex}/{strain}_{demultiplex}_filtered.fastq'
	conda:
		'.....yml'
	params:
		outputdir = '{path}/{main_folder}/preprocessing/{demultiplex}/',
		prefix = '{strain}_{demultiplex}_filtered'
	threads: 8
	shell:
		''


# using nanoplot for quality statistics of the long reads
rule nanoplot:
	input:
		rules.filtlong.output.reads
	output:
		'{path}/{main_folder}/quality/nanoplot/{strain}/{strain}_{demultiplex}NanoPlot-report.html'
	conda:
		'read_quality.yml'
	params:
		outputdir = '{path}/{main_folder}/quality/nanoplot/{strain}/',
		prefix = '{strain}_{demultiplex}'
	threads: 8
	shell:
		'NanoPlot -t {threads} --minlength 1000 --fastq {input} -o {params.outputdir} -p {params.prefix} --title {params.prefix}_minlength1000'

# long read assembler
# Flye
rule flye:
	input:
		'{path}/{main_folder}/preprocessing/{demultiplex}/{strain}_{demultiplex}.fastq'
	output:
		contigs = '{path}/{main_folder}/assembly/{strain}_{demultiplex}_flye/assembly.fasta'
	conda:
		'flye.yml'
	params:
		outputdir = '{path}/{main_folder}/assembly/{strain}_{demultiplex}_flye/',
	threads: 32
	shell:
		'flye --asm-coverage 50 --nano-raw {input} -o {params.outputdir} -t {threads} -g 1000000'
        # --asm-coverage 50         coverage: 50X
        # -g 1000000                genome size 1 million

# rename flye output
rule rename_flye:
	input:
		contigs = rules.flye.output.contigs
	output:
		flye = '{path}/{main_folder}/assembly/{strain}_{demultiplex}_flye/{strain}_{demultiplex}_flye.fasta'
	shell:
		'mv {input.contigs} {output.flye}'

# polishing the assembly with Racon four times using the long reads		use minimap2 for the mapping inbetween the polishing
rule minimap2_racon_long:
	input:
		assembly = '{path}/{main_folder}/assembly/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}.fasta',
		reads = '{path}/{main_folder}/preprocessing/{demultiplex}/{strain}_{demultiplex}.fastq'
	output:
		out = '{path}_docker/racon/{main_folder}/{strain}_{demultiplex}_{assembler}_long4.fasta'
	params:
		strain = '{strain}',
		demultiplex = '{demultiplex}',
		assembler = '{assembler}',
		path = '{path}/{main_folder}',
		main_folder = '{main_folder}'
	threads: 32
	run:
		import os
		import time
		reads_path = str(params.path) + '/preprocessing/' + str(params.demultiplex) + '/'
		only_reads = str(params.strain) + '_' + str(params.demultiplex) + '.fastq'
        # change the docker path please
		racon_path = '/mnt/prostlocal2/projects/st_mycoplasma_assembly_docker/racon/' + str(params.main_folder) + '/'
		assembly_path = str(params.path) + '/assembly/' + str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '/'
		paf_path = str(params.path) + '/postprocessing/' + str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '/'
		assembly_file = str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '.fasta'
		for i in range(4):
			if i == 0:
				minimap2_input_assembly = assembly_path + assembly_file
				racon_in_assembly = assembly_file
			else:
				minimap2_input_assembly = racon_path + out_assembly
				racon_in_assembly = out_assembly
			out_assembly = str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '_long' + str(i + 1) + '.fasta'
			paf = str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '_long' + str(i + 1) + '.paf'
			minimap2 = 'minimap2 -x map-ont -t 16 ' + minimap2_input_assembly + ' ' + reads_path + only_reads + ' > ' + paf_path + paf
			print(minimap2 + '\n')
			os.system(minimap2)
			while os.path.isfile(paf_path + paf) == False:
				time.sleep(5)
			if i == 0:
				racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + paf_path + ':/input2 -v ' + assembly_path + ':/input3 -v ' + racon_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t 32 /input1/' + only_reads + ' /input2/' + paf + ' /input3/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
			else:
				racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + paf_path + ':/input2 -v ' + racon_path + ':/input3 -v ' + racon_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t 32 /input1/' + only_reads + ' /input2/' + paf + ' /input3/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
			print(racon + '\n')			
			os.system(racon)
		
rule move_racon_long:
	input:
		racon_out = rules.minimap2_racon_long.output.out
	output:
		move = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_long4.fasta'
	shell:
		'mv {input.racon_out} {output.move}'

# polishing the assembly with medaka
rule medaka:
	input:
		reads = '{path}/{main_folder}/preprocessing/{demultiplex}/{strain}_{demultiplex}.fastq',
		racon_out = rules.move_racon_long.output.move
	output:
		medaka = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/consensus.fasta'
	conda:
		'postprocessing.yml'
	params:
		outputdir = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/'
	threads: 32
	shell:
		'medaka_consensus -i {input.reads} -d {input.racon_out} -o {params.outputdir} -t {threads} -m r941_min_high'

# polishing the assembly with Racon four times using the short reads		use minimap2 for the mapping inbetween the polishing
rule minimap2_racon_short:
	input:
		assembly = rules.medaka.output.medaka,
		reads = rules.unique_readID.output.unique
	output:
		out = '{path}_docker/racon/{main_folder}/{strain}_{demultiplex}_{assembler}_short4.fasta'
	params:
		strain = '{strain}',
		demultiplex = '{demultiplex}',
		assembler = '{assembler}',
		path = '{path}/{main_folder}',
		main_folder = '{main_folder}'
	threads: 32
	run:
		import os
		import time
		reads_path = str(params.path) + '/preprocessing/'
		only_reads = str(params.strain) + '_unique.fastq'
        # change the docker path please
		racon_path = '/mnt/prostlocal2/projects/st_mycoplasma_assembly_docker/racon/' + str(params.main_folder) + '/'
		assembly_path = str(params.path) + '/postprocessing/' +  str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '/'
		paf_path = str(params.path) + '/postprocessing/' + str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '/'
		assembly_file = 'consensus.fasta'
		for i in range(4):
			if i == 0:
				minimap2_input_assembly = assembly_path + assembly_file
				racon_in_assembly = assembly_file
			else:
				minimap2_input_assembly = racon_path + out_assembly
				racon_in_assembly = out_assembly
			out_assembly = str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '_short' + str(i + 1) + '.fasta'
			paf = str(params.strain) + '_' + str(params.demultiplex) + '_' + str(params.assembler) + '_short' + str(i + 1) + '.paf'
			minimap2 = 'minimap2 -x sr -t 16 ' + minimap2_input_assembly + ' ' + reads_path + only_reads + ' > ' + paf_path + paf
			print(minimap2 + '\n')
			os.system(minimap2)
			while os.path.isfile(paf_path + paf) == False:
				time.sleep(5)
			if i == 0:
				racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + paf_path + ':/input2 -v ' + assembly_path + ':/input3 -v ' + racon_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t 32 /input1/' + only_reads + ' /input2/' + paf + ' /input3/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
			else:
				racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + paf_path + ':/input2 -v ' + racon_path + ':/input3 -v ' + racon_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t 32 /input1/' + only_reads + ' /input2/' + paf + ' /input3/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
			print(racon + '\n')			
			os.system(racon)

# move Racon output from docker folder to work folder
rule move_racon_short:
	input:
		racon_out = rules.minimap2_racon_short.output.out
	output:
		move = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_short4.fasta'
	shell:
		'mv {input.racon_out} {output.move}'

rule circlator:
	input:
        assembly = rules.move_racon_short.output.move
    output:
        circ = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}_circ.fasta'
    conda:
        '....yml'
    # params:
    #     outputdir = '{path}/{main_folder}/postprocessing/{strain}_{demultiplex}_{assembler}/',
    #     prefix = '{strain}_{demultiplex}_{assembler}_circ'
    threads: 16
    shell:
        ''

rule prokka:
    input:
        assembly = rules.circlator.output.circ
    output:
        gff = '{path}/{main_folder}/prokka/{strain}_{demultiplex}_{assembler}/{strain}_{demultiplex}_{assembler}.gff'
    conda:
        'prokka.yml'
    params:
        outputdir = '{path}/{main_folder}/prokka/{strain}_{demultiplex}_{assembler}/',
        prefix = '{strain}_{demultiplex}_{assembler}'
    threads: 16
    shell:
        'prokka --cpus {threads} --gcode 4 --outdir {params.outputdir} --force --prefix {params.prefix} {input.assembly}'

rule quast:
	input:
        # insert path to final assembly file to incorporate all in quast statistics
		in1 = '{path}/{main_folder}/postprocessing/STRAIN1_qcat_flye/STRAIN1_qcat_flye_short4.contigs.fasta',
		in2 = '{path}/{main_folder}/postprocessing/STRAIN2_qcat_flye/STRAIN2_qcat_flye_short4.contigs.fasta',
		in3 = '{path}/{main_folder}/postprocessing/STRAIN3_qcat_flye/STRAIN3_qcat_flye_short4.contigs.fasta'
		# in1 = '{path}/{main_folder}/postprocessing/14DD0147_qcat_flye/14DD0147_qcat_flye_short4.contigs.fasta',
		# in2 = '{path}/{main_folder}/postprocessing/15DD0207_qcat_flye/15DD0207_qcat_flye_short4.contigs.fasta',
		# in3 = '{path}/{main_folder}/postprocessing/15DD0165_qcat_flye/15DD0165_qcat_flye_short4.contigs.fasta',
		# in4 = '{path}/{main_folder}/postprocessing/15DD0234_qcat_flye/15DD0234_qcat_flye_short4.contigs.fasta',
		# in5 = '{path}/{main_folder}/postprocessing/15DD0163_qcat_flye/15DD0163_qcat_flye_short4.contigs.fasta',
		# in6 = '{path}/{main_folder}/postprocessing/15DD0233_qcat_flye/15DD0233_qcat_flye_short4.contigs.fasta',
		# in7 = '{path}/{main_folder}/postprocessing/15DD0164_qcat_flye/15DD0164_qcat_flye_short4.contigs.fasta',
		# in8 = '{path}/{main_folder}/postprocessing/15DD0218_qcat_flye/15DD0218_qcat_flye_short4.contigs.fasta',
		# in9 = '{path}/{main_folder}/postprocessing/15DD0210_qcat_flye/15DD0210_qcat_flye_short4.contigs.fasta',
		# in10 = '{path}/{main_folder}/postprocessing/15DD0238_qcat_flye/15DD0238_qcat_flye_short4.contigs.fasta',
		# in11 = '{path}/{main_folder}/postprocessing/15DD0261_qcat_flye/15DD0261_qcat_flye_short4.contigs.fasta',
		# in12 = '{path}/{main_folder}/postprocessing/15DD0228_qcat_flye/15DD0228_qcat_flye_short4.contigs.fasta'
	output:
		report = '{path}/{main_folder}/quality/quast/report.html'
		#...
	conda:
		'assembly_quality.yml'
	params:
		outputdir = '{path}/{main_folder}/quality/quast/',
		path = '{path}/mycoplasma_bovis_genomes'
	threads: 16
	shell:
		'quast {input.in1} {input.in2} {input.in3} -r {params.path}/mycoplasma_bovis_referenceGenome.fasta -g {params.path}/GCF_000183385.1_ASM18338v1_genomic.gff -o {params.outputdir}'
		# {input.in4} {input.in5} {input.in6} {input.in7} {input.in8} {input.in9} {input.in10}