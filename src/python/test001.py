#!/usr/bin/python
import numpy as np
import scipy.sparse as spm
import pylab as plt
import scipy.sparse.linalg as spla
from scipy import sparse 


 
def fun1(K,Kreg,R,Rl):
    normK = np.max(np.abs(K.diagonal()))
    n1 = np.linalg.norm(sparse.csc_matrix.dot(K,R.toarray()))#/normK
    n2 = np.linalg.norm(sparse.csc_matrix.dot(K.transpose(),Rl.toarray()))#/normK
#    #
    print('||K  * R1||=',n1)
    print('||Kt * R2||=',n2)
    print('||K||=',normK)
#    
    #K_full=K.toarray()
    
    iKreg = spla.splu(Kreg)
    MP_cond = np.zeros((K.shape[0],K.shape[0]))
    for i in range(K.shape[0]):
        MP_cond[:,i] = \
                sparse.csc_matrix.dot(K,iKreg.solve(K[:,i].toarray()[:,0]))\
                -K[:,i].toarray()[:,0]
                      
    print("||K-K*pinv(K)*K|| =",np.abs(MP_cond).sum())
    return 0
 
 
for i in range(len(mat_K)):
    for j in range(len(mat_K[i])):
        out_ = fun1(mat_K[i][j],mat_Kreg[i][j],mat_R1[i][j],mat_R2[i][j])
        
        
        
        
