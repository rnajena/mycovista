#create folders for the output files from the nanopore data pipeline
import os

configfile = open('config.yaml', 'r').read()
liste = configfile.split('\n')
if liste[0][:5] == 'path:':
    path = liste[0][7:len(liste[0]) - 1] + '/hybrid'
    os.mkdir(path)
    liste.pop(0)
else:
    path = os.getcwd()
# os.mkdir(path + '/raw_data')
for i in range (len(liste)):
    parameter = liste[i].split(':')
    details = parameter[1]
    parameter = parameter[0]
    details = details.split(',')
    details[0] = details[0][2:len(details[0])]
    details[len(details) - 1] = details[len(details) - 1][0:len(details[len(details) - 1]) - 1]
    # if parameter != 'strains':
    #     os.mkdir(path + '/' + parameter)
    if parameter == 'strains':
        strains = details
    # if parameter == 'demultiplexing':
    #     demultiplexing = details
    if parameter == 'assembly':
        os.mkdir(path + '/' + parameter)
        assembly = details
    # if parameter == 'polishing':
    #     polishing = details
    if parameter == 'quality':
        os.mkdir(path + '/' + parameter)
        quality = details
# for preprocesser in demultiplexing:
#     os.mkdir(path + '/demultiplexing' + '/' + preprocesser)
for bacteria in strains:
    for assembler in assembly:
        os.mkdir(path + '/assembly' + '/' + bacteria + '_' + assembler)
            # for polisher in polishing:
            #     os.mkdir(path + '/polishing' + '/' + bacteria + '_' + preprocesser + '_' + assembler + '_' + polisher)
for check in quality:
    os.mkdir(path + '/quality' + '/' + check)