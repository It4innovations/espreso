
import numpy as np 
from scipy import sparse 
import myModul as mM
import config_espreso_python  
import pylab as plt
import scipy.sparse.linalg as spla

n_clus          = 2
n_subPerClust   = 4



problem_info = {'n_clus': n_clus,'n_subPerClust':n_subPerClust}

path = '../log/'

mat_K       = []
mat_Salfa   = []
mat_Kreg    = []
mat_B0      = []
mat_B1      = []
mat_R       = []
vec_f       = []
vec_c       = []
vec_weight  = []
vec_index_weight = []
mat_Schur   = []
for i in range(n_clus): 
    mat_K.append([])
    mat_Kreg.append([])
    mat_B0.append([])
    mat_B1.append([])
    mat_R.append([])
    vec_f.append([])
    vec_c.append([])
    vec_weight.append([])
    vec_index_weight.append([])
    mat_Salfa.append(mM.load_matrix(path,'Salfa',0,'',makeSparse=False,makeSymmetric=True))
    for j in range(n_subPerClust):  
        mat_K[i].append(mM.load_matrix(path,'K',i,j,makeSparse=True,makeSymmetric=False))
        mat_Kreg[i].append(mM.load_matrix(path,'Kreg',i,j,makeSparse=True,makeSymmetric=True))      
        mat_B0[i].append(mM.load_matrix(path,'B0',i,j,makeSparse=True,makeSymmetric=False))
        mat_B1[i].append(mM.load_matrix(path,'B1',i,j,makeSparse=True,makeSymmetric=False))
        mat_R[i].append(mM.load_matrix(path,'R',i,j,makeSparse=True,makeSymmetric=False))
        vec_f[i].append(mM.load_vector(path,'f',i,j))
        vec_c[i].append(mM.load_vector(path,'c',i,j))
        vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        tmp = mM.load_vector(path,'loc_ind_weight',i,j)
        tmp = tmp.astype(np.int32)
        vec_index_weight[i].append(tmp)
        
###############################################################################
####################### FETI PREPROCESSING ####################################
###############################################################################   



conf = config_espreso_python

if conf.precondDualSystem=='dirichlet':
    for i in range(len(mat_K)):
        mat_Schur.append([])
        for j in range(len(mat_K[i])):
            I = np.arange(0,mat_K[i][j].shape[0])
            _B1 = mat_B1[i][j].copy()
            J = np.unique(_B1.tocoo().col)            
            I = np.delete(I,J)
            K_II = mat_K[i][j][np.array(I)[:,np.newaxis],np.array(I)].tocsc()           
            K_IJ = mat_K[i][j][np.array(I)[:,np.newaxis],np.array(J)].toarray()
            K_JJ = mat_K[i][j][np.array(J)[:,np.newaxis],np.array(J)].toarray()
            iK_II = spla.splu(K_II)
            iK_II_K_IJ = np.zeros(K_IJ.shape)
            for k in range(K_IJ.shape[1]):
                iK_II_K_IJ[:,k] = iK_II.solve(K_IJ[:,k])
            mat_Schur[i].append([J,K_JJ-np.dot(K_IJ.transpose(),iK_II_K_IJ)])


if True:
    mat_B0ker   = []
    
    i_B0ker     = []
    j_B0ker     = []
    v_B0ker     = []
    for i in range(len(mat_B0)):

        mat_B0ker.append([])
        i_B0ker.append([])
        j_B0ker.append([])
        v_B0ker.append([])
        
        for j in range(len(mat_B0[i])-1):
            for k in range(j+1,len(mat_B0[i])):
                tmpB_jk = sparse.hstack([np.abs(mat_B0[i][j]),np.abs(mat_B0[i][k])])
                indx = np.ravel(tmpB_jk.sum(1)==2) 

                if (np.sum(indx)>0):
                                        
                    iB0_j_bool = np.zeros(mat_B0[i][j].shape[1],dtype=bool)
                    ij         = mat_B0[i][j].tocsr().indices
                    
                    iB0_k_bool = np.zeros(mat_B0[i][k].shape[1],dtype=bool)
                    ik         = mat_B0[i][k].tocsr().indices
                    
                    iB0_j = mat_B0[i][j].tocsr()[indx,:].indices
                    iB0_k = mat_B0[i][k].tocsr()[indx,:].indices

                    R_g_j = mat_R[i][j].toarray()[iB0_j,:]
                    R_g_k = mat_R[i][k].toarray()[iB0_k,:]
                    
                    print(R_g_j-R_g_k)

                    

#plt.clf()

u,lam = mM.feti(mat_K,mat_Kreg,vec_f,mat_Schur,mat_B1,vec_c,vec_weight,\
                            vec_index_weight,mat_R)
                            
                            
uHDP,lamH = mM.hfeti(mat_K,mat_Kreg,vec_f,mat_Schur,mat_B0,mat_B1,vec_c,\
                        vec_weight,vec_index_weight,mat_R,mat_Salfa)
#
norm_del_u = 0
norm_u = 0
for i in range(len(u)):
    for j in range(len(u[i])):
        norm_del_u += np.linalg.norm(u[i][j]-uHDP[i][j])
        norm_u += np.linalg.norm(u[i][j])

print('|u_TFETI-u_HTFETI|/|u_TFETI| = ',norm_del_u/norm_u)


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
