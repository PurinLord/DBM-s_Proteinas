import numpy as np
import gzip
import cPickle
import sys
import random
from getAttrfromSeq import *
import seqVectorizer as sv

class dataProssesor(object):

	#def __init__(self):
		##inicio
	
	def loadCathDomainSeqs(self, data):
		
		print "cargando " + data

		#se abre el archivo con las etiquetas y las secuencias de aminoacidos
		#fseqs = open("CathDomainSeqs.ATOM.v3.5.0")
		fseqs = open(data)
		seqs = []
		labels = []
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

		return seqs, labels


	def load26VectorRepresentation(self, data):

		#se carga el archivo que tiene los vectores de 26 y las etiquetas de las secuencias
		print "cargando " + data
		dic26 = {}
		with gzip.open(data, 'rn') as f:
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

		return dic26


	def joinData(self, dataSet, dictionary):

		#ahora se unen la entrada y la salida
		print "Relacionando la entrada con la salida"
		protA26 = []
		#vector con las entradas que no tengan una imagen 3D correspondeinte
		#para el preprosesamiento
		preProsses = []
		#print len(dataSet[1])
		for i in range(len(dataSet[1])):
			try:
				#usando la etiquta se busca dada elemento de dataSet en el diciionario de 26
				#si lo encuentra lo une a la salida
				protA26.append([dataSet[0][i], dictionary[dataSet[1][i]]])
			except ValueError:
				#si no se encuentra se guarda para el preprossesamiento
				preProsses.append([dataSet[0][i], 0])
				#print "NOT"
		#print len(protA26)
		#print len(preProsses)
		#print protA26
		return protA26, preProsses


	def splitData(self, mainData, preProssesData = [], ratio = 7):

		random.shuffle(mainData)
		random.shuffle(preProssesData)
		#se calcula una septima parte de la longutdud
		#para tomar el vector de entreneamiento el de prieba y el de validacion
		subSet = len(mainData)/ratio
		valid = mainData[0:subSet]
		test = mainData[subSet:subSet * 2]
		train = mainData[subSet * 2: -1] 
		return train, valid, test, preProssesData

	def formating (self, data):
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
	
	def packing (self, data, fileName):
		f = gzip.open(fileName, 'w')
		f.write(cPickle.dumps(data, 1))
		
		f.close()

	def standardFormatAndPack(self, train, valid, test, preProsses = [], verbos = False):
		print len(train)
		vecTrain = self.formating(train)
		vecVal = self.formating(valid)
		vecTest = self.formating(test)
		#vector con las cadenas sin correspondencia de 26 para el pre entrenamiento
		vectFullTrain = self.formating(train + preProsses)
		
		if verbos == True:
			print len(vecTrain)
			print len(vecTrain[0])
			print type(vecTrain[0][0][0])
			np.set_printoptions(threshold=np.nan)
			print "len entrada " + str(len(vecTrain[0][0][0]))
			print vecTrain[0][0][0]
			print "len salida " + str(len(vecTrain[0][1][0]))
			print vecTrain[0][1][0]
		
		for i in range(len(vecTrain)):
			self.packing([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")
		
		#empacando la lista del pre prosesamiento
		#se le agregan vectores 0 en el lugar de valid y test para que DBN.py las pueda usar directamente
		a = []
		for i in range(2):
			a.append([0,0])
		self.packing([vectFullTrain[0],a ,a], "Pre.pkl.gz")


	def standardProsses (self, data1, data2):
		seqs, labels = self.loadCathDomainSeqs(data1)
		dic26 = self.load26VectorRepresentation(data2)

		vectorizer = sv.seqVectorizer(comblength=7)
		vectNormal = vectorizer.fit_transform(seqs)
		protSeq = [vectNormal, labels]

		mainData, preProsses = self.joinData(protSeq, dic26)
		
		train, valid, test, preProsses = self.splitData(mainData, preProsses) 
		print len(train)
		print len(preProsses)

		self.standardFormatAndPack(train, valid, test, preProsses, True)
def main():
	print "comienzo"
	_dataProssesor = dataProssesor()
	_dataProssesor.standardProsses("CathReducida", "CATHALL.txt.tar.gz")

if __name__== '__main__':
    main()
