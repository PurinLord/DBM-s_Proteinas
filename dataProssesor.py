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
                #print len(sallable) ,
        
                #crea el vector principal de 26 dimenciones
                data26 = map(int, vector[1:27])
        
                #print sallable
                #print data26
        
                #se cea el diccionario 
                dic26[sallable] = data26

        return dic26


    def joinData(self, dataSet, dictionary):
        """
        Matches a data set with amino acid sequence with its corresponding 26 dimensional vector
        Returns: matrix with matching sequence-26vector, vector all the seq that didn't have a patch

<<<<<<< HEAD
		#ahora se unen la entrada y la salida
		print "Relacionando la entrada con la salida"
		print type(dictionary)

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
				preProsses.append([dataSet[0][i], 0])
				#print "NOT"
		print cont
		print "Prot a 26"
		print len(protA26)
		print "pre train"
		print len(preProsses)

		return protA26, preProsses
=======

        :type dataSet: list 
        :param dataSet: list of lables and aminoacid seguence
        :type dictionary: dict
        :param dictionary: dict of  lables and 26 vector
        """
>>>>>>> 01cd8e03241395e7b4e7b43ad8da953e6cf17993

        #ahora se unen la entrada y la salida
        print "Relacionando la entrada con la salida"
        print type(dictionary)
        protA26 = []
        #vector con las entradas que no tengan una imagen 3D correspondeinte
        #para el preprosesamiento
        preProsses = []
        #print len(dataSet[1])
        for i in range(len(dataSet[1])):
            D = dictionary.get(dataSet[1][i], -1)

            if D > 0:
                #usando la etiquta se busca dada elemento de dataSet en el diciionario de 26
                #si lo encuentra lo une a la salida
                protA26.append([dataSet[0][i], dictionary[dataSet[1][i]]])
            else:
                #si no se encuentra se guarda para el preprossesamiento
                preProsses.append([dataSet[0][i], 0])
                #print "NOT"
        #print len(protA26)
        #print len(preProsses)
        #print protA26
        print
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
        random.shuffle(preProssesData)
        #se calcula una septima parte de la longutdud
        #para tomar el vector de entreneamiento el de prieba y el de validacion
        subSet = len(mainData)/ratio
        valid = mainData[0:subSet]
        test = mainData[subSet:subSet * 2]
        train = mainData[subSet * 2: -1] 
        return train, valid, test, preProssesData

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

    def standardFormatAndPack(self, train, valid, test, preProsses = [], verbos = False):
        """
        Converts the 3 data sets (train, calid, test) in to 26 data sets in de format udes by DBM.py
        and one for the pre prosesing
        """
        print len(train)
        print len(train[0])
        print len(train[0][0])
        vecTrain = self.formating(train)
        vecVal = self.formating(valid)
        vecTest = self.formating(test)
        #vector con las cadenas sin correspondencia de 26 para el pre entrenamiento
        vectFullTrain = self.formating([train[0] + preProsses, train[1] + preProsses])
        
        if verbos == True:
            print len(vecTrain)
            print len(vecTrain[0])
            print type(vecTrain[0][0][0])
            np.set_printoptions(threshold=np.nan)
            print "len entrada " + str(len(vecTrain[0][0][0]))
            print vecTrain[0][0][0]

<<<<<<< HEAD
		random.shuffle(mainData)
		random.shuffle(preProssesData)
		#se calcula una septima parte de la longutdud
		#para tomar el vector de entreneamiento el de prieba y el de validacion
		subSet = len(mainData)/ratio
		valid = mainData[0:subSet]
		test = mainData[subSet:subSet * 2]
		train = mainData[subSet * 2: -1] 
		print "train, val, test, prePr"
		print len(train), len(valid), len(test), len(preProssesData)
		return train, valid, test, preProssesData
	"""
	Puts the data in the correcto format for DBM.py
	"""
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

	def formatingOne (self, data):
		out = []
		#print len(data)
		#print len(data[0][1])
		#print len(data[0][0])
		#se revisa que data[0][1] no sea un entero, para el caso de preProsses
		print type(data[0][1][0])
		entrada = []
		for i in range(len(data)):
			entrada.append(data[i][0])
			salida = data[i][1] #= 0
		out = [entrada, salida]
		return out
	"""
	Packs the data in the format used by DBM.py using the chones name
	"""
	def packing (self, data, fileName):
		f = gzip.open(fileName, 'w')
		f.write(cPickle.dumps(data, 1))
		
		f.close()
	"""
	Converts the 3 data sets (train, calid, test) in to 26 data sets in de format udes by DBM.py
	and one for the pre prosesing
	"""
	def standardFormatAndPack(self, train, valid, test, preProsses = [], verbos = False):

		vecTrain = self.formating(train)
		vecVal = self.formating(valid)
		vecTest = self.formating(test)
		#vector con las cadenas sin correspondencia de 26 para el pre entrenamiento
		vectFullTrain = self.formatingOne(train + preProsses)

		if verbos == True:
			print "trian"
			print len(train), len(train[0]), len(train[0][0]), len(train[0][1])
			print "peProsses"
			print len(preProsses), len(preProsses[0]), len(preProsses[0][0])

			print "vTrain"
			print len(vecTrain), len(vecTrain[0]), len(vecTrain[0][0]), len(vecTrain[0][0][0]) 
			print "vPre"
			print len(vectFullTrain), len(vectFullTrain[0]), len(vectFullTrain[0][0])
		
			print type(vecTrain[0][0][0])
			np.set_printoptions(threshold=np.nan)
		
		for i in range(len(vecTrain)):
			self.packing([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")
		
		#empacando la lista del pre prosesamiento
		#se le agregan vectores 0 en el lugar de valid y test para que DBN.py las pueda usar directamente
		a = []
		for i in range(2):
			a.append([0,0])
		self.packing([vectFullTrain[0],a ,a], "Pre.pkl.gz")
=======
            print vecTrain[0][1][0]
        
        for i in range(len(vecTrain)):
            self.packing([vecTrain[i], vecVal[i], vecTest[i]], "ProteinData-" + str(i + 1) + ".pkl.gz")
        
        #empacando la lista del pre prosesamiento
        #se le agregan vectores 0 en el lugar de valid y test para que DBN.py las pueda usar directamente
        a = []
        for i in range(2):
            a.append([0,0])
        self.packing([vectFullTrain[0],a ,a], "Pre.pkl.gz")

>>>>>>> 01cd8e03241395e7b4e7b43ad8da953e6cf17993

    def standardProsses (self, data1, data2, comblength):
        """
        Prosseses de file "data1" and the file "data2" using all the methods and converting the aminoacid
        sequence with seqVectorizer
        Returns: 26 data sest to be used un DBM.py and one for the pre-training
        """
        seqs, labels = self.loadCathDomainSeqs(data1)
        dic26 = self.load26VectorRepresentation(data2)

        vectorizer = sv.seqVectorizer(comblength=comblength)
        print str(len(seqs)) + " transformando secuencias"
        #print seqs
        vectNormal = vectorizer.fit_transform(seqs)
        
        #Packing of seqVectorizer object. This file is NEEDED to transform new sequences
        with open('vectorizer.pkl', 'wb') as fp:
            cPickle.dump(vectorizer, fp)
        
        protSeq = [vectNormal, labels]

<<<<<<< HEAD
		mainData, preProsses = self.joinData(protSeq, dic26)
		
		train, valid, test, preProsses = self.splitData(mainData, preProsses) 
=======
        mainData, preProsses = self.joinData(protSeq, dic26)
        
        train, valid, test, preProsses = self.splitData(mainData, preProsses) 
        print len(train)
        print "prep" 
        print len(preProsses)
>>>>>>> 01cd8e03241395e7b4e7b43ad8da953e6cf17993

        self.standardFormatAndPack(train, valid, test, preProsses, True)
def main():
    print "comienzo"
    _dataProssesor = dataProssesor()
    #_dataProssesor.standardProsses("CathReducida", "CATHALL.txt.tar.gz")
    _dataProssesor.standardProsses("CathDomainSeqs.ATOM.v3.5.0", "CATHALL.txt.tar.gz", 2)

if __name__== '__main__':
    main()
