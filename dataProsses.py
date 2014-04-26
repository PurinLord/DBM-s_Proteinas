import numpy as np
import gzip
import cPickle
import sys
import random
from getAttrfromSeq import *
import seqVectorizer as sv

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

salSeq = []
vect26 = []
vectClasi = []

#Lee el archivo en gz
print "cargando CATHALL.txt.tar.gz"
with gzip.open("CATHALL.txt.tar.gz", 'rn') as f:
	#lee cada linea del archivo
	for line in f:
		#print line
		#transforma la linea en un vector
		vector = line.split(',')
		#elimina el ultimo caracter (\n)
		vector[-1] = vector[-1][:-1]
		#print vector
		#crea el vector principal de 26 dimenciones
		data26 = map(int, vector[1:27])
		#print data26
		#crea el resto con las etiquetas
		sallable = [vector[0]]
		#print dominioYclas
		#se va creando el vector con todos los vectorsitos de 26
		vect26.append(data26)
		#lo mismo para las etiquetas
		vectClasi.append(sallable[0])

#####datos de finales de salida
salSeq = [vect26[:-1], vectClasi[:-1]]
#se elimina el ultimo elemento porque esta vacio, no se porque
#print salSeq[0][1]
#print salSeq[1][1]

#ahora se unen la entrada y la salida
print "Relacionando la entrada con la salida"
protA26 = []
print len(protSeq[1])
for i in range(len(protSeq[1])):
	try:
		#dada la estiqueta de la entrada se busca su correspondiente en la salida
		index = salSeq[1].index(protSeq[1][i])
		#se guarda en protA26 la entrada con su salida correspondiente
		protA26.append([protSeq[0][i], salSeq[0][index]])
	except ValueError:
		print "NOT"
print len(protA26)
#print protA26

random.shuffle(protA26)
#se calcula una septima parte de la longutdud
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

for i in range(len(vecTrain)):
	empacar([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")


	

