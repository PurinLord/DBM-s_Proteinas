from getAttrfromSeq import *
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.ndimage import convolve
from sklearn import linear_model, datasets, metrics
from sklearn.cross_validation import train_test_split
from sklearn.neural_network import BernoulliRBM
from sklearn.pipeline import Pipeline
from sklearn import linear_model
from sklearn import svm

normalAAslist = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
F_Ic4list = ['AWM','GST','HPY','CVIFL','DNQ','ER','K']
LESKlist = ['AST', 'CVILWYMPF','HQN','RK','ED','G']
MSlist = ['AVLIMC','WYHF','TQSN','RK','ED','GP']
Lmap =  F_Ic4list

import cPickle
RBMfilename = sys.argv[3]#'RBMtrained.txt'
rbm = cPickle.load(open(RBMfilename, 'rb'))
print 'trained RBM loaded'

#plt.imshow(rbm.components_)
#plt.savefig('RBMtmp.png',dpi=450)
#plt.show()

fvectors = open(sys.argv[2])
Dvectors = dict()
for v in fvectors:
	v = v.strip().split(',')
	domname = v[0]
	vector = map(int,v[1:-1])
	Dvectors[domname] = vector


fseqs = open(sys.argv[1])
transformedSeqs = []
domname = '_'
DtransformedSeq = dict()
comblength = 7
for s in fseqs:
	s = s.strip()
	if s[0] == '>':
		domname = s[8:15] 
		continue
	if 'X' in s: continue
	if not domname in Dvectors.keys(): continue
	transformedSeq = rbm.transform(np.array(createAAFreqVector(s,Lmap,comblength)))
	DtransformedSeq[domname] = transformedSeq
	#transformedSeqs.append(transformedSeq)
	
#Now there are dictionaries with transformed sequences and vectors associated with each domain
#Lets put this in a matrix

vectors = []
transformedSeqs = []
for k in Dvectors.iterkeys():
	if not k in DtransformedSeq: continue #this should not happen, but...
	vectors.append(Dvectors[k])
	transformedSeqs.append(DtransformedSeq[k])
	
vectors = np.asarray(vectors)	
transformedSeqs = np.asarray(transformedSeqs)
print 'vectors shape' , vectors.shape
print 'transformedSeqs shape' , transformedSeqs.shape
	
#ok, now it is time to train the regressor
#lets make a list of regressors

from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score	

numdims = 26 #Number of regressors, just for testing
regresorList = []
for i in xrange(numdims):
	#clf = linear_model.SGDRegressor(loss="squared_epsilon_insensitive",n_iter=1000,eta0=0.005)
	#clf = clf = linear_model.Ridge(alpha=0.005)
	#clf = linear_model.LinearRegression()
	#clf = linear_model.BayesianRidge()
	#clf = linear_model.Lasso()
	#clf = linear_model.ARDRegression()
	#clf = linear_model.OrthogonalMatchingPursuit()
	#clf = linear_model.LassoLars()
	#clf = linear_model.ElasticNet()
	#clf = linear_model.Perceptron()
	#clf = linear_model.PassiveAggressiveClassifier()
	#clf  = svm.SVR()
	clf = svm.SVR(kernel='rbf', C=1e5, gamma=0.1) #C=1e5 funciona muy nice
	clf.fit(transformedSeqs, vectors[:,i])
	regresorList.append(clf) #Tfor save memory
	#ytrue = vectors[:,i]
	#ypredict = np.asarray(map(int,clf.predict(transformedSeqs)))
	#for x in zip(ytrue,ypredict):
	#	print x
	#print explained_variance_score(ytrue, ypredict)	, r2_score(ytrue, ypredict)
	

RegresorsListFN = 'RegresorsList264SVMNormallist.txt'
print 'Initiating dump of list of regressors'
cPickle.dump(regresorList, open(RegresorsListFN, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)
print 'dumping process finished'

#for i in xrange(numdims):
#	ytrue = vectors[:,i]
#	ypredict = np.asarray(map(int,regresorList[i].predict(transformedSeqs)))
#	print ytrue.shape, ypredict.shape
#	print explained_variance_score(ytrue, ypredict)	
	#for j in xrange(len(transformedSeqs)):
		#print regresorList[i].predict(transformedSeqs[j]), vectors[j][i] 
	
	
#plt.imshow(transformedSeqs[:150])
#plt.savefig('transformedSeqs.png',dpi=450)
#plt.show()


