import itertools
from collections import Counter
import numpy as np

def dictFromListMapping(L):
	"""
	creates a dictionary containing the mapping depicted in list L
	"""
	D = dict()
	for i in xrange(len(L)):
		for a in L[i]:
			D[a] = str(i)
	return D



def mapSeq(S,M):
	"""
	Returns the transformed string S given the mapping in dictionary M
	"""
	res = ''.join([M[c] for c in S])
	return res


def getFreqsDict(seq,dictmap,n):
	"""
	returns a list with frequences of combinations from 2 to n of aminoacids in seq with residue mapping dictmap
	"""
	combsCounter = Counter()
	for i in xrange(2,n+1):
		for j in xrange(len(seq)-i):
			strKey = ''.join(sorted(mapSeq(seq[j:j+i],dictmap)))
			combsCounter[strKey] += 1
			
	return combsCounter
	
def iterAACombs(n,alfabet):
	"""
	returns a list with all combinations of aminoacids from length 2 to n
	"""	
	#AAs = 'ARNDCEQGHILKMFPSTWYV'
	#F_Ic4 = '1234567'
	#AAs = '1234567'
	AAs = alfabet
	#print AAs
	AAcombsList = []
	for i in xrange(2,n+1):
		for combs in itertools.combinations_with_replacement(AAs,i): #itertools.product(AAs, repeat=i): 
			yield ''.join(sorted(combs))


def createAAFreqVector(seq,L,n,normalize=False):
	"""
	returns an array with frequecies of each aa combination
	"""
	import string
	d = dictFromListMapping(L)
	AAfreqsList = []
	C = getFreqsDict(seq,d,n)
	alfab = ''.join([string.ascii_letters[i] for i in xrange(len(L))])
	for i in iterAACombs(n,alfab):
		AAfreqsList.append(C[i])
	if normalize: return  (AAfreqsList - np.min(AAfreqsList, 0)) / (np.max(AAfreqsList, 0) + 0.0001)
	return AAfreqsList

	
#a = 'MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQD'
#b = getFreqsDict(a,4)
#print b

#AAFL = createAAFreqVector(a,4)
#print len(AAFL)
#for i in AAFL:
#	if i:
#		print i
