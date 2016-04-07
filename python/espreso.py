
import numpy as np 
from scipy import sparse 
import myModul as mM
import config_espreso_python  

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
#vec_weight  = []
for i in range(n_clus): 
    mat_K.append([])
    mat_Kreg.append([])
    mat_B0.append([])
    mat_B1.append([])
    mat_R.append([])
    vec_f.append([])
    mat_Salfa.append(mM.load_matrix(path,'Salfa',0,'',makeSparse=False,makeSymmetric=True))
    #vec_weight.append([])
    for j in range(n_subPerClust):  
        mat_K[i].append(mM.load_matrix(path,'K',i,j,makeSparse=True,makeSymmetric=False))
        mat_Kreg[i].append(mM.load_matrix(path,'Kreg',i,j,makeSparse=True,makeSymmetric=True))      
        mat_B0[i].append(mM.load_matrix(path,'B0',i,j,makeSparse=True,makeSymmetric=False))
        mat_B1[i].append(mM.load_matrix(path,'B1',i,j,makeSparse=True,makeSymmetric=False))
        mat_R[i].append(mM.load_matrix(path,'R',i,j,makeSparse=True,makeSymmetric=False))
        vec_f[i].append(mM.load_vector(path,'f',i,j))
        #vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        
 
for i in range(n_clus):
    for j in range(n_subPerClust):
        if (i==0 and j==0): 
            f       = vec_f[0][0]
            #weight  = vec_weight[0][0] 

        else:            
            f       = np.concatenate((f,vec_f[i][j]))        
            #weight  = np.concatenate((weight,vec_weight[i][j]))          
                      



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

weight = 1
u,lam = mM.feti(mat_K,mat_Kreg,f,mat_B1,mat_R,weight)

uHDP,lamH = mM.hfeti(mat_K,mat_Kreg,f,B0,B1,R,mat_Salfa,weight)



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
