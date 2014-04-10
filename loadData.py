import gzip
import cPickle
import readline

archive = "CATHALL.txt.tar.gz"

salida = []
vect26 = []
vectclasi = []

#Lee el archivo en gz
with gzip.open(archive, 'rn') as f:
	#lee cada linea del archivo
	for line in f:
		#print line
		#transforma la linea en un vector
		vector = line.split(',')
		#elimina el ultimo caracter (\n)
		vector[-1] = vector[-1][:-1]
		#print vector
		#crea el vector principal de 26 dimenciones
		data26 = vector[1:27]
		#print data26
		#crea el resto con las etiquetas
		dominioYclasi = [vector[0], vector[-1]]
		#print dominioYclas
		#se va creando el vector con todos los vectorsitos de 26
		vect26.append(data26)
		#lo mismo para las etiquetas
		vectclasi.append(dominioYclasi)

#se crea el vector con ambas cosas
#como el que se le alimenta a la DBM
train= [vect26, vectclasi]
#las otras dos salidas se copian de la primera para que tenga el mismo formato que las databeses de logistic_sgd.py
valid = [vect26, vectclasi]
test = [vect26, vectclasi]
#se crea el vector de salida
salida = [salida, valid, test]

#se comprime con cPickle y luego con gzip para que tenga el mismo formato que usa la DBM
f = gzip.open('datasetProt.pkl.gz', 'w')
f.write(cPickle.dumps(salida, 1))

f.close()


print len(salida)
print len(salida[0]) 
