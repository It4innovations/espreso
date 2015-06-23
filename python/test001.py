import numpy as np
import pylab as plt

nS  = 10

K = []; f = []; B0 = []
print('... reading data ...')
for i in range(nS):
    Kloc = np.loadtxt('../dumped_files/K_mat_'+str(i))
    K.append(Kloc)
    floc = np.loadtxt('../dumped_files/f_vec_'+str(i))
    f.append(floc)
    B0loc = np.loadtxt('../dumped_files/B0_vec_'+str(i))
    B0.append(B0loc)

