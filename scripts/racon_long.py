import os
import time

reads_path = snakemake.params[3] + '/preprocessing/' + snakemake.params[1] + '/'
only_reads = snakemake.input[1].split('/')[len(snakemake.input[1].split('/'))-1]
assembly_path = snakemake.params[3] + '/assembly/' + snakemake.params[0] + '_' + snakemake.params[1] + '_' + snakemake.params[2] + '/'
out_path = snakemake.params[3] + '/postprocessing/' + snakemake.params[0] + '_' + snakemake.params[1] + '_' + snakemake.params[2] + '/'
assembly_file = snakemake.input[0].split('/')[len(snakemake.input[0].split('/'))-1]
for i in range(4):
	if i == 0:
		minimap2_input_assembly = assembly_path + assembly_file
		racon_in_assembly = assembly_file
	else:
		minimap2_input_assembly = out_path + out_assembly
		racon_in_assembly = out_assembly
	out_assembly = snakemake.params[0] + '_' + snakemake.params[1] + '_' + snakemake.params[2] + '_long' + str(i + 1) + '.fasta'
	paf = snakemake.params[0] + '_' + snakemake.params[1] + '_' + snakemake.params[2] + '_long' + str(i + 1) + '.paf'
	minimap2 = 'minimap2 -x map-ont -t ' + str(snakemake.threads) + ' ' + minimap2_input_assembly + ' ' + reads_path + only_reads + ' > ' + out_path + paf
	print(minimap2 + '\n')
	os.system(minimap2)
	while os.path.isfile(out_path + paf) == False:
		time.sleep(5)
	if i == 0:
		racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + out_path + ':/input2 -v ' + assembly_path + ':/input3 -v ' + out_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t ' + str(snakemake.threads) + ' /input1/' + only_reads + ' /input2/' + paf + ' /input3/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
	else:
		racon = 'docker run --rm --user $(id -u):$(id -g) -it -v ' + reads_path + ':/input1 -v ' + out_path + ':/input2 -v ' + out_path + ':/output quay.io/biocontainers/racon:1.3.2--he941832_0 sh -c "racon -t ' + str(snakemake.threads) + ' /input1/' + only_reads + ' /input2/' + paf + ' /input2/' + racon_in_assembly + ' > /output/' + out_assembly + '"'
	print(racon + '\n')			
	os.system(racon)
