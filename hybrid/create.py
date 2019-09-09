#create folders for the output files from the nanopore data pipeline
import os

configfile = open('config.yaml', 'r').read()
liste = configfile.split('\n')
for i in liste:
    if i == '' or i[0] == '#':
        liste.remove(i)
if liste[0][:5] == 'path:':
    path = liste[0][7:len(liste[0]) - 1] + '/hybrid'
    os.mkdir(path)
    liste.pop(0)
else:
    path = os.getcwd()
os.mkdir(path + '/raw_data')
for i in range(len(liste)):
    parameter = liste[i].split(':')
    details = parameter[1]
    parameter = parameter[0]
    details = details.split(',')
    details[0] = details[0][2:len(details[0])]
    details[len(details) - 1] = details[len(details) - 1][0:len(details[len(details) - 1]) - 1]
    if parameter == 'strains':
        strains = details
    if parameter == 'preprocessing_short':
        os.mkdir(path + '/preprocessing')
        preprocessing_short = details
    if parameter == 'demultiplexing':
        demultiplexing = details
    if parameter == 'assembly':
        os.mkdir(path + '/' + parameter)
        assembly = details
    if parameter == 'postprocessing':
        os.mkdir(path + '/' + parameter)
        postprocessing = details
    if parameter == 'quality':
        os.mkdir(path + '/' + parameter)
        quality = details
for preprocesser in preprocessing_short:
    os.mkdir(path + '/preprocessing' + '/' + preprocesser)
for demultiplexer in demultiplexing:
    os.mkdir(path + '/preprocessing' + '/' + demultiplexer)
    for bacteria in strains:
        for assembler in assembly:
            os.mkdir(path + '/assembly' + '/' + bacteria + '_' + demultiplexer + '_' + assembler)
            os.mkdir(path + '/postprocessing' + '/' + bacteria + '_' + demultiplexer + '_' + assembler)
for check in quality:
    os.mkdir(path + '/quality' + '/' + check)
    for bacteria in strains:
        os.mkdir(path + '/quality' + '/' + check + '/' + bacteria)