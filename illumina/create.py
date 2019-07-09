#create Illumina folders for the output files
import os

configfile = open('config.yaml', 'r').read()
liste = configfile.split('\n')
if liste[0][:5] == 'path:':
    path = liste[0][7:len(liste[0]) - 1]
    os.mkdir(path)
    liste.pop(0)
else:
    path = os.getcwd()
for i in range (len(liste)):
    parameter = liste[i].split(':')
    details = parameter[1]
    parameter = parameter[0]
    details = details.split(',')
    help = details[0]
    details[0] = help[2:len(help)]
    help = details[len(details) - 1]
    details[len(details) - 1] = help[0:len(help) - 1]
    if parameter != 'strains':
        os.mkdir(path + '/' + parameter)
    if parameter == 'strains':
        strains = details
    if parameter == 'preprocessing':
        preprocessing = details
    if parameter == 'assembly':
        assembly = details
    if parameter == 'quality':
        quality = details
os.mkdir(path + '/raw_data')
os.mkdir(path + '/fastqc')
os.mkdir(path + '/fastqc/preprocessed')
for bacterium in strains:
    os.mkdir(path + '/preprocessing/' + bacterium)
    for assembler in assembly:
        os.mkdir(path + '/assembly/' + bacterium + '_' + assembler)
        for check in quality:
            os.mkdir(path + '/quality/' + check)