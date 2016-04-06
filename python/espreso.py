
import numpy as np 
from scipy import sparse 
import myModul as mM
import config_espreso_python  

n_clus          = 1
n_subPerClust   = 2

path = '../log/'

mat_K       = []
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
    #vec_weight.append([])
    for j in range(n_subPerClust):  
        mat_K[i].append(mM.load_matrix(path,'K',i,j,makeSparse=True,makeSymmetric=False))
        mat_Kreg[i].append(mM.load_matrix(path,'Kreg',i,j,makeSparse=True,makeSymmetric=True))      
        mat_B0[i].append(mM.load_matrix(path,'B0',i,j,makeSparse=True,makeSymmetric=False))
        mat_B1[i].append(mM.load_matrix(path,'B1',i,j,makeSparse=True,makeSymmetric=False))
        mat_R[i].append(mM.load_matrix(path,'R',i,j,makeSparse=False,makeSymmetric=False))
        vec_f[i].append(mM.load_vector(path,'f',i,j))
        #vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        
 
for i in range(n_clus):
    for j in range(n_subPerClust):
        if (i==0 and j==0):
            B0      = mat_B0[0][0].copy()  
            B1      = mat_B1[0][0].copy()      
            K       = mat_K[0][0]   
            Kreg    = mat_Kreg[0][0]
            R       = mat_R[0][0]
            f       = vec_f[0][0]
            #weight  = vec_weight[0][0]
            diagR   = np.sum(mat_R[0][0]*mat_R[0][0],axis=1)

        else: 
            if (B0.__class__.__name__ != 'list'):
                B0      = sparse.hstack((B0,mat_B0[i][j]))
            B1      = sparse.hstack((B1,mat_B1[i][j]))
            K       = sparse.block_diag((K,mat_K[i][j]))
            Kreg    = sparse.block_diag((Kreg,mat_Kreg[i][j]))
            R       = sparse.block_diag((R,mat_R[i][j]))            
            f       = np.concatenate((f,vec_f[i][j]))        
            #weight  = np.concatenate((weight,vec_weight[i][j]))
            diagR   = np.concatenate((diagR,np.sum(mat_R[i][j]*mat_R[i][j],axis=1)))           
            
diagR   = diagR 
if n_clus*n_subPerClust==1:
    R       = sparse.csc_matrix(R) 
weight = 1
            
mat_K       = []
mat_Kreg    = []
mat_B0      = []
mat_B1      = []
mat_R       = []
vec_f       = []
vec_weight  = []            



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


u,lam = mM.feti(K,Kreg,f,B1,R,diagR,weight)

uHDP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,diagR,weight)



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
