import os
import time

reads_path = snakemake.params[1] + '/preprocessing/long_reads/'
only_reads = snakemake.input[1].split('/')[len(snakemake.input[1].split('/'))-1]
assembly_path = snakemake.params[1] + '/assembly/' + snakemake.params[0] + '_flye/'
out_path = snakemake.params[1] + '/postprocessing/' + snakemake.params[0] + '_flye/'
assembly_file = snakemake.input[0].split('/')[len(snakemake.input[0].split('/'))-1]

for i in range(4):
	if i == 0:
		minimap2_input_assembly = assembly_path + assembly_file
		racon_in_assembly = assembly_file
	else:
		minimap2_input_assembly = out_path + out_assembly
		racon_in_assembly = out_assembly

	out_assembly = snakemake.params[0] + '_flye_long' + str(i + 1) + '.fasta'
	paf = snakemake.params[0] + '_flye_long' + str(i + 1) + '.paf'
	minimap2 = 'minimap2 -x map-ont -t ' + str(snakemake.threads) + ' ' + minimap2_input_assembly + ' ' + reads_path + only_reads + ' > ' + out_path + paf
	
	print(minimap2 + '\n')
	os.system(minimap2)

	while os.path.isfile(out_path + paf) == False:
		time.sleep(5)

	if i == 0:
		racon = 'racon -t ' + str(snakemake.threads) + ' ' + reads_path + only_reads + ' ' + out_path + paf + ' ' + assembly_path + racon_in_assembly + ' > ' + out_path + out_assembly
	else:
		racon = 'racon -t ' + str(snakemake.threads) + ' ' + reads_path + only_reads + ' ' + out_path + paf + ' ' + out_path + racon_in_assembly + ' > ' + out_path + out_assembly
	
	print('\n' + racon + '\n')			
	os.system(racon)
