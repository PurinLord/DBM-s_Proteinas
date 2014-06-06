import numpy as np
import gzip
import cPickle
import sys
import random
from getAttrfromSeq import *
import seqVectorizer as sv


class dataProssesor(object):
    """    This class can read specific data files and unpak  
    its content it can process data vectors format and pack them into the the type used by dbm.py
    it also contains a methos thoes all the data formating in one step

    loadCathDomainSeqs(data)
    load26VectorRepresentation(data)
    joinData(dataSet, dictionary):
    splitData(mainData, preProssesData, ratio):
    formating (data):
    packing (data, fileName):
    standardFormatAndPack(train, valid, test, preProsses, verbos):
    standardProsses (data1, data2):
    
    Reads a file of plain text containing data in the format:
    >domain|"id lable"|3_5_0
    "amino acid sequence" (with no separation)
    Example
    >domain|1ahqA00|3_5_0
    GIAVSDDCVQKFNELKLGHQHRYVTFKMNASNTEVVVEHVGGPNATYEDFKSQLPERDCRYAIFDYEFQVDGGQRNKITFILWAPDSAPIKSKMMYTSTKDSIKKKLVGIQVEVQATDAAEISEDAVSERAKK
    Returns: vector containing sequence, vector containing labels

    :type data: string
    :param data: name of the file to unpak
    """
    def loadCathDomainSeqs(self, data):
        
        print "cargando " + data

        #se abre el archivo con las etiquetas y las secuencias de aminoacidos
        fseqs = open(data)
        seqs = []
        labels = []
        for s in fseqs:
            
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
        """ 
        Reads a file of the type .txt.tar.gz containing data in the format:
        "id lable","26 dimention vector"(separated by colons),#_#_##_#
        Example
        1ivsA04,2,0,2,6,6,6,25,4,0,0,2,1,0,6,3,0,0,0,0,0,0,0,0,0,0,5,3_30_1170_10
        Returns: dictionary with lables and vectors

        :type data: string
        :param data: name of the file to unpak
        """

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
        
                #crea el vector principal de 26 dimenciones
                data26 = map(int, vector[1:27])
        
                #se cea el diccionario 
                dic26[sallable] = data26

        return dic26


    def joinData(self, dataSet, dictionary):
        """
        Matches a data set with amino acid sequence with its corresponding 26 dimensional vector
        Returns: matrix with matching sequence-26vector, vector all the seq that didn't have a patch


        :type dataSet: list 
        :param dataSet: list of lables and aminoacid seguence
        :type dictionary: dict
        :param dictionary: dict of  lables and 26 vector
        """

        #ahora se unen la entrada y la salida

	print "data"
	print len(dataSet) , len(dataSet[0])
	print "dic"
	print len(dictionary)

        protA26 = []
        #vector con las entradas que no tengan una imagen 3D correspondeinte
        #para el preprosesamiento
        preProsses = []
        #print len(dataSet[1])
	cont = 0
        for i in range(len(dataSet[1])):
            D = dictionary.get(dataSet[1][i], -1)

            if D > 0:
                #usando la etiquta se busca dada elemento de dataSet en el diciionario de 26
                #si lo encuentra lo une a la salida
                protA26.append([dataSet[0][i], dictionary[dataSet[1][i]]])
            else:
                #si no se encuentra se guarda para el preprossesamiento
		cont += 1 
                preProsses.append([dataSet[0][i], [0]*26])
		
	print "no match"
	print cont
	print "Prot a 26"
	print len(protA26)
	print "pre train"
	print len(preProsses)

        return protA26, preProsses

    def splitData(self, mainData, preProssesData = [], ratio = 7):
        """
        Split the data set in to the training, validation an test set
        it also creates a data set for the pre-training useing preProssesData
        Returns: training, validation, test sets, and teh pre-training set
          
        :type mainData: list 
        :param mainData: list of aminoacid sequence and 26 vector 
        :type preProssesData
        :param preProssesData
        :type ratio
        :param ratio
        """

        random.shuffle(mainData)
        #se calcula una septima parte de la longutdud
        #para tomar el vector de entreneamiento el de prieba y el de validacion
        subSet = len(mainData)/ratio
        valid = mainData[0:subSet]
        test = mainData[subSet:subSet * 2]
        train = mainData[subSet * 2: -1] 
	fullPreProsses = train + preProssesData 
        random.shuffle(fullPreProsses)
	print "train, val, test, prePr"
	print len(train), len(valid), len(test), len(fullPreProsses)
        return train, valid, test, fullPreProsses

    def formating (self, data):
        """
        Puts the data in the correcto format for DBM.py
        """
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
        """
        Packs the data in the format used by DBM.py using the chones name
        """
        f = gzip.open(fileName, 'w')
        f.write(cPickle.dumps(data, 1))
        
        f.close()

	def packAll(data1, data2, comblength):
		
		seqs, labels = self.loadCathDomainSeqs(data1)
		dic26 = self.load26VectorRepresentation(data2)
		vectorizer = sv.seqVectorizer(comblength=comblength)
		vectNormal = vectorizer.fit_transform(seqs)
		protSeq = [vectNormal, labels]

		mainData, preProsses = self.joinData(protSeq, dic26)
		print len(mainData)
		print len(mainData[0])
		print len(mainData[0][0])
		print len(mainData[0][1])

		print len(preProsses)
		print len(preProsses[0])
		print len(preProsses[0][0])
		print len(preProsses[0][1])
		
		train, valid, test, preProsses = self.splitData(mainData, preProsses) 
		print len(train)
		print len(train[0])
		print len(train[0][0])
		print len(train[0][1])

		print len(preProsses)
		print len(preProsses[0])
		print len(preProsses[0][0])
		print len(preProsses[0][1])

		vecTrain = self.formating(train)
		vecVal = self.formating(valid)
		vecTest = self.formating(test)

		vectFullTrain = self.formating(preProsses)
		
		print "trian"
		print len(train), len(train[0]), len(train[0][0]), len(train[0][1])
		print "peProsses"
		print len(preProsses), len(preProsses[0]), len(preProsses[0][0])

		print "vTrain"
		print len(vecTrain), len(vecTrain[0]), len(vecTrain[0][0]), len(vecTrain[0][0][0]) 
		print "vPre"
		print len(vectFullTrain), len(vectFullTrain[0]), len(vectFullTrain[0][0]), len(vectFullTrain[0][0][0])

		for i in range(len(vecTrain)):
		    self.packing([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")
		
		#empacando la lista del pre prosesamiento
		#se le agregan vectores 0 en el lugar de valid y test para que DBN.py las pueda usar directamente
		a = []
		for i in range(2):
		    a.append([0,0])
		self.packing([vectFullTrain[0], vecVal[0], vecTest[0]], "Pre.pkl.gz")

def main():
	print "comienzo"
	_dataProssesor = dataProssesor()
	data1 = "CathDomainSeqs.ATOM.v3.5.0"
	data2 = "CATHALL.txt.tar.gz"
	comblength = 2
	
        seqs, labels = _dataProssesor.loadCathDomainSeqs(data1)
        dic26 = _dataProssesor.load26VectorRepresentation(data2)
        vectorizer = sv.seqVectorizer(comblength=comblength)
        vectNormal = vectorizer.fit_transform(seqs)
        protSeq = [vectNormal, labels]

        mainData, preProsses = _dataProssesor.joinData(protSeq, dic26)
	print len(mainData)
	print len(mainData[0])
	print len(mainData[0][0])
	print len(mainData[0][1])

	print len(preProsses)
	print len(preProsses[0])
	print len(preProsses[0][0])
	print len(preProsses[0][1])
        
        train, valid, test, preProsses = _dataProssesor.splitData(mainData, preProsses) 
	print len(train)
	print len(train[0])
	print len(train[0][0])
	print len(train[0][1])

	print len(preProsses)
	print len(preProsses[0])
	print len(preProsses[0][0])
	print len(preProsses[0][1])

        vecTrain = _dataProssesor.formating(train)
        vecVal = _dataProssesor.formating(valid)
        vecTest = _dataProssesor.formating(test)

	vectFullTrain = _dataProssesor.formating(preProsses)
        
	print "trian"
	print len(train), len(train[0]), len(train[0][0]), len(train[0][1])
	print "peProsses"
	print len(preProsses), len(preProsses[0]), len(preProsses[0][0])

	print "vTrain"
	print len(vecTrain), len(vecTrain[0]), len(vecTrain[0][0]), len(vecTrain[0][0][0]) 
	print "vPre"
	print len(vectFullTrain), len(vectFullTrain[0]), len(vectFullTrain[0][0]), len(vectFullTrain[0][0][0])

        for i in range(len(vecTrain)):
            _dataProssesor.packing([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")
        
        #empacando la lista del pre prosesamiento
        #se le agregan vectores 0 en el lugar de valid y test para que DBN.py las pueda usar directamente
        a = []
        for i in range(2):
            a.append([0,0])
        _dataProssesor.packing([vectFullTrain[0], vecVal[0], vecTest[0]], "Pre.pkl.gz")

if __name__== '__main__':
	main()

