import numpy as np
import gzip
import cPickle
import sys
import random
from getAttrfromSeq import *
import seqVectorizer as sv

#se abre el archivo con las etiquetas y las secuencias de aminoacidos
#fseqs = open("CathDomainSeqs.ATOM.v3.5.0")
fseqs = open("CathReducida")
seqs = []
labels = []
#se lee cad alinea del archivo de secuencias
print "cargando CathDomainSeqs.ATOM.v3.5.0"
for s in fseqs:
	
	#print s
	s = s.strip()

	#si comienza con ">" es la etiqueta de una secuencia
	#entonces se guarda en el vector secuencia
	if s[0] == '>':
		split = s.split('|')

		labels.append(split[1])
	else:
	#de lo contrario es la secuencia que sigue a la etiqueta
		if 'X' in s: 
	#si la secuencia contienen "X" no es valida, asi que se borra su etiqueta correspondiente
			del labels[-1]
			continue 
		seqs.append(s)

#print len(seqs)

vectorizer = sv.seqVectorizer(comblength=7)

vectNormal = vectorizer.fit_transform(seqs)


######datos finales de entrada
protSeq = [vectNormal, labels]
dicSeq = dict(zip(labels, vectNormal))


dic26 = {}
#se carga el archivo que tiene los vectores de 26 y las etiquetas de las secuencias
print "cargando CATHALL.txt.tar.gz"
with gzip.open("CATHALL.txt.tar.gz", 'rn') as f:
	#lee cada linea del archivo
	for line in f:

		#transforma la linea en un vector
		vector = line.split(',')
		#elimina el ultimo caracter (\n)
		vector[-1] = vector[-1][:-1]

		#crea el vector con las etiquetas
		sallable = vector[0]
		#la primera entrada no e sun dato
		#por seguridad si la etiqueta es mas larga que 7 se descarta la entrada
		if len(sallable) > 7:
			continue
		#print len(sallable) ,

		#crea el vector principal de 26 dimenciones
		data26 = map(int, vector[1:27])

		#print sallable
		#print data26

		#se cea el diccionario 
		dic26[sallable] = data26

#####datos de finales de salida

#ahora se unen la entrada y la salida
print "Relacionando la entrada con la salida"
protA26 = []
#vector con las entradas que no tengan una imagen 3D correspondeinte
#para el preprosesamiento
preProsses = []
print len(protSeq[1])
for i in range(len(protSeq[1])):
	try:
		#usando la etiquta se busca dada elemento de protSeq en el diciionario de 26
		#si lo encuentra lo une a la salida
		protA26.append([protSeq[0][i], dic26[protSeq[1][i]]])
	except ValueError:
		#si no se encuentra se guarda para el preprossesamiento
		preProsses.append([protSeq[0][i], 0])
		print "NOT"
print len(protA26)
print len(preProsses)
#print protA26

random.shuffle(protA26)
random.shuffle(preProsses)
#se calcula una septima parte de la longutdud
#para tomar el vector de entreneamiento el de prieba y el de validacion
subSet = len(protA26)/7
valid = protA26[0:subSet]
test = protA26[subSet:subSet * 2]
train = protA26[subSet * 2: -1] 

#print "train"
#print len(train)
#print len(train[0])
#print len(train[0][0])
#print len(train[0][1])
#print type(train[0][1])
#print train[0][1]
#print train[0][0]


def format (data):
	out = []
	#print len(data)
	#print len(data[0][1])
	#print len(data[0][0])
	for j in range(len(data[0][1])):
		entrada = []
		salida = []
		for i in range(len(data)):
			entrada.append(data[i][0])
			salida.append(data[i][1][j])
		out.append([entrada, salida])
	return out

def empacar (data, fileName):
	f = gzip.open(fileName, 'w')
	f.write(cPickle.dumps(data, 1))
	
	f.close()

#prueba = format(train)
#print "prueba"
#print len(prueba)
#print len(prueba[0])
#print len(prueba[0][0])
#print len(prueba[0][1])
#print len(prueba[0][0][0])
#print prueba[0][1]
#print len(prueba[1][0][0])

vecTrain = format(train)
vecVal = format(valid)
vecTest = format(test)
#vector con las cadenas sin correspondencia de 26 para el pre entrenamiento
vectFullTrain = format(train + preProsses)

for i in range(len(vecTrain)):
	empacar([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")

#empacando la lista del pre prosesamiento
empacar([vectFullTrain, 0, 0], "Pre.pkl.gz")

	

