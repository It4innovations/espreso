
import numpy as np 
from scipy import sparse 
import myModul as mM
import config_espreso_python  
import pylab as plt

n_clus          = 2
n_subPerClust   = 2



problem_info = {'n_clus': n_clus,'n_subPerClust':n_subPerClust}

path = '../log/'

mat_K       = []
mat_Salfa   = []
mat_Kreg    = []
mat_B0      = []
mat_B1      = []
mat_R       = []
vec_f       = []
vec_weight  = []
vec_index_weight = []
#vec_weight  = []
for i in range(n_clus): 
    mat_K.append([])
    mat_Kreg.append([])
    mat_B0.append([])
    mat_B1.append([])
    mat_R.append([])
    vec_f.append([])
    vec_weight.append([])
    vec_index_weight.append([])
    mat_Salfa.append(mM.load_matrix(path,'Salfa',0,'',makeSparse=False,makeSymmetric=True))
    #vec_weight.append([])
    for j in range(n_subPerClust):  
        mat_K[i].append(mM.load_matrix(path,'K',i,j,makeSparse=True,makeSymmetric=False))
        mat_Kreg[i].append(mM.load_matrix(path,'Kreg',i,j,makeSparse=True,makeSymmetric=True))      
        mat_B0[i].append(mM.load_matrix(path,'B0',i,j,makeSparse=True,makeSymmetric=False))
        mat_B1[i].append(mM.load_matrix(path,'B1',i,j,makeSparse=True,makeSymmetric=False))
        mat_R[i].append(mM.load_matrix(path,'R',i,j,makeSparse=True,makeSymmetric=False))
        vec_f[i].append(mM.load_vector(path,'f',i,j))
        vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        tmp = mM.load_vector(path,'loc_ind_weight',i,j)
        tmp = tmp.astype(np.int32)
        vec_index_weight[i].append(tmp)
        #vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        
###############################################################################
####################### FETI PREPROCESSING ####################################
###############################################################################            
#np.hstack;np.vstack;np.column_stack;np.row_stack


#FETI=False

#if FETI:
    #B = sparse.vstack((B0 ,B1 ))
    #u,lam = mM.feti(K,Kreg,f,B1,R,weight)
#else:    
    #u,lam = mM.hfeti(K,Kreg,f,B0,B1,R,weight)

#B = sparse.vstack((B0 ,B1 ))
#u,lam = mM.feti(K,Kreg,f,B,R,weight)

conf = config_espreso_python


plt.clf()

u,lam = mM.feti(mat_K,mat_Kreg,vec_f,mat_B1,vec_weight,\
                            vec_index_weight,mat_R)
                            
                            
uHDP,lamH = mM.hfeti(mat_K,mat_Kreg,vec_f,mat_B0,mat_B1,vec_weight,\
                        vec_index_weight,mat_R,mat_Salfa)
#
#norm_del_u = 0
#norm_u = 0
#for i in range(len(u)):
#    for j in range(len(u[i])):
#        norm_del_u += np.linalg.norm(u[i][j]-uHDP[i][j])
#        norm_u += np.linalg.norm(u[i][j])
#
#print('|u_TFETI-u_HTFETI|/|u_TFETI| = ',norm_del_u/norm_u)


#conf.iterative_Kplus=False
#uHDP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,weight)
#
#methodToImproveSolByIterMethod      = 'cg_dx'
#conf.precondPrimalSystem            = 'diag' 
#conf.iterative_Kplus=True
#uHSP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,weight)



#ndu = np.linalg.norm(uHSP-uHDP)
#nu = np.linalg.norm(uHDP)
#print('||uHDP-uHSP||/||uHDP||       = ',ndu/nu)

#ndu = np.linalg.norm(u-uH)
#nu = np.linalg.norm(u)
#ndlam = np.linalg.norm(lam-lamH)
#nlam = np.linalg.norm(lam)
#print('||u-uH||/||u||       = ',ndu/nu)
#print('||lam-lamH||/||lam|| = ',ndlam/nlam)
