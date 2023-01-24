# create folders for output_folder files generated during the de novo hybrid assembly pipeline

import os
import glob

configfile = open('config.yaml', 'r').read()
liste = configfile.split('\n')

for elem in liste[:]:
    if elem == '' or elem[0] == '#':
        liste.remove(elem)
        continue
    config = elem.split(': [')[1][:-1]
    if 'pipeline:' in elem: pipeline = config
    elif 'input_shortreads:' in elem: input_short = config
    elif 'input_longreads:' in elem: input_long = config
    elif 'output:' in elem:
        output_folder = config
        if os.path.exists(output_folder) == False:
            os.mkdir(output_folder)
        if len(output_folder) == 0:
            output_folder = os.getcwd()
    elif 'strains:' in elem:
        strains = elem.split(': ')[1][1:-1].split(',')

os.mkdir(f"{output_folder}/annotation")

os.mkdir(f"{output_folder}/assembly")

os.mkdir(f"{output_folder}/final_assemblies")

os.mkdir(f"{output_folder}/preprocessing")
os.mkdir(f"{output_folder}/preprocessing/long_reads")
os.mkdir(f"{output_folder}/postprocessing")

os.mkdir(f"{output_folder}/quality")

os.mkdir(f"{output_folder}/quality/nanoplot")
os.mkdir(f"{output_folder}/quality/quast")

os.mkdir(f"{output_folder}/raw_data")

if pipeline == 'hybrid':
    os.mkdir(f"{output_folder}/preprocessing/short_reads")
    os.mkdir(f"{output_folder}/preprocessing/short_reads/fastp")
    os.mkdir(f"{output_folder}/preprocessing/short_reads/trimmomatic")
    os.mkdir(f"{output_folder}/quality/fastqc")

for elem in strains:
    os.mkdir(f"{output_folder}/annotation/{elem}_flye")
    os.mkdir(f"{output_folder}/assembly/{elem}_flye")
    os.mkdir(f"{output_folder}/postprocessing/{elem}_flye")
    os.mkdir(f"{output_folder}/quality/nanoplot/{elem}")

    if pipeline == 'hybrid':
        for file in glob.glob(f"{input_short}/*{elem}*"):
            if '_1.' in file: os.system(f"ln -s {file} {output_folder}/raw_data/{elem}_1.fastq.gz")
            elif '_2.' in file: os.system(f"ln -s {file} {output_folder}/raw_data/{elem}_2.fastq.gz")
    for file in glob.glob(f"{input_long}/*{elem}*"):
        if '_1.' not in file and '_2.' not in file:
            os.system(f"ln -s {file} {output_folder}/raw_data/{elem}.fastq")



   