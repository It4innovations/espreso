#!/usr/bin/python
import numpy as np
import scipy.sparse as spm
import pylab as plt
import scipy.sparse.linalg as sla

nS  = 10

K = []; f = []; B0 = []
print('... reading data ...')
for i in range(nS):
    tmp0 = np.loadtxt('../dumped_files/K_mat_'+str(i))
    tmp1 = spm.coo_matrix((tmp0[1::,2],(np.int32(tmp0[1::,0]-1),np.int32(tmp0[1::,1]-1))), shape=(tmp0[0,0],tmp0[0,1])).tocsr()
    K.append(tmp1)
    tmp0 = np.loadtxt('../dumped_files/B0_mat_'+str(i))
    tmp1 = spm.coo_matrix((tmp0[1::,2],(np.int32(tmp0[1::,0]-1),np.int32(tmp0[1::,1]-1))), shape=(tmp0[0,0],tmp0[0,1])).tocsr()
    B0.append(tmp1)
    tmp0 = np.loadtxt('../dumped_files/f_vec_'+str(i))
    f.append(tmp0)
    print('%d-th sub.' % i)

P = []; L = []; U = []
print(' lu decomposition')
for i in range(nS):
    Kplus_Bt = spm.linalg.splu(K[i],B0[0].transpose())
