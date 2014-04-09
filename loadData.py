import gzip
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
salida= [vect26, vectclasi]
print len(salida)
print len(salida[0]) 
