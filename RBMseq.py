from getAttrfromSeq import *
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.ndimage import convolve
from sklearn import linear_model, datasets, metrics
from sklearn.cross_validation import train_test_split
from sklearn.neural_network import BernoulliRBM
from sklearn.pipeline import Pipeline

normalAAslist = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
F_Ic4list = ['AWM','GST','HPY','CVIFL','DNQ','ER','K']
MSlist = ['AVLIMC','WYHF','TQSN','RK','ED','GP']
LESKlist = ['AST', 'CVILWYMPF','HQN','RK','ED','G']
Lmap =  F_Ic4list #MSlist #LESKlist #F_Ic4list
#MS = dictFromListMapping(MSlist)
#F_Ic4 = dictFromListMapping(F_Ic4list)
#print F_Ic4
#normalAAs = dictFromListMapping(normalAAslist)
###############################################################################
# Setting up




# Models we will use
logistic = linear_model.LogisticRegression()
rbm = BernoulliRBM(random_state=0, verbose=True,batch_size=15)

classifier = Pipeline(steps=[('rbm', rbm), ('logistic', logistic)])

###############################################################################
# Training

# Hyper-parameters. These were set by cross-validation,
# using a GridSearchCV. Here we are not performing cross-validation to
# save time.
rbm.learning_rate = 0.001
rbm.n_iter = 50
# More components tend to give better prediction performance, but larger
# fitting time
rbm.n_components = 250
#logistic.C = 6000.0

fseqs = open(sys.argv[1])
seqs = []
for s in fseqs:
	s = s.strip()
	if s[0] == '>': continue
	if 'X' in s: continue
	seqs.append(s)
#print seqs

#seqs = ['MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQD',
#'QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAEKMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTSVLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHPFLFLIKHNPTNTIVYFGRYWSP',
#'MVDREQLVQKARLAEQAERYDDMAAAMKNVTELNEPLSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTSADGNEKKIEMVRAYREKIEKELEAVCQDVLSLLDNYLIKNCSETQYESKVFYLKMKGDYYRYLAEVATGEKRATVVESSEKAYSEAHEISKEHMQPTHPIRLGLALNYSVFYYEIQNAPEQACHLAKTAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDD',
#'MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSL',
#'MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD',
#]

comblength = 7

X = map(lambda s : np.array(createAAFreqVector(s,Lmap,comblength)) , seqs)
#print X

#X = (X - np.min(X, 0)) / (np.max(X, 0) + 0.0001)  # 0-1 scaling
#print X.shape

rbm.fit(X)
ssss ='MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSL'
transformedSeq = rbm.transform(np.array(createAAFreqVector(ssss,Lmap,comblength)))
print transformedSeq
print 'len', len(transformedSeq)
# Training RBM-Logistic Pipeline
#classifier.fit(X_train, Y_train)

# Training Logistic regression
#logistic_classifier = linear_model.LogisticRegression(C=100.0)
#logistic_classifier.fit(X_train, Y_train)

###############################################################################
# Evaluation

print()

###############################################################################
# Plotting
#for i in (rbm.components_):
#	print len(i)

#for i, comp in enumerate(rbm.components_):
#	print i, len(comp)
print rbm.components_.shape

#plt.imshow(rbm.components_)
#plt.savefig('RBMtmp.png',dpi=450)
#plt.figure(figsize=(4.2, 4))
#for i, comp in enumerate(rbm.components_):
#    plt.subplot(10, 10, i + 1)
#    plt.imshow(comp, cmap=plt.cm.gray_r,
#               interpolation='nearest')
#    plt.xticks(())
#    plt.yticks(())
#plt.suptitle('100 components extracted by RBM', fontsize=16)
#plt.subplots_adjust(0.08, 0.02, 0.92, 0.85, 0.08, 0.23)

plt.show()


#Dumping RBM
import cPickle
RBMfilename = 'RBMtrained'+str(rbm.n_components)+'neuronsFIc4list.txt'
print 'Initiating dumps'
cPickle.dump(rbm, open(RBMfilename, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL) 
print 'Dumping done'
