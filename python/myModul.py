import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import config_espreso_python as conf
import pylab as plt
#import sys  
#
#
def load_matrix(path,str0,i,j,makeSparse,makeSymmetric): 
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile, ndmin=2)
    if (tmp.shape[0]==1):
        tmp = []
    else:
        n = np.int32(tmp[0,0])   
        m = np.int32(tmp[0,1])
        I = tmp[1::,0]-1;    J = tmp[1::,1]-1;    V = tmp[1::,2]
#    
        print(str0,i,j)
        if (makeSymmetric):
            logInd = J>I; 
            I = np.concatenate((I,J[logInd]))
            J = np.concatenate((J,I[logInd]))
            V = np.concatenate((V,V[logInd]))    
        
        if (makeSparse):
            tmp = sparse.csc_matrix((V,(I,J)),shape=(n,m)).tocsc()
        else:
            if (m==1):
                tmp = V
            else:                
                tmp = sparse.csc_matrix((V,(I,J)),shape=(n,m)).toarray()
#
    return tmp
###############################################################################
def load_vector(path,str0,i,j):
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile)
    return tmp 
###############################################################################
class DENS_SOLVE: 
    # solving of A*x = b with A = L*Lt
    # x = Lt^-1 * ( L^-1 * b)
    def __init__(self, A):
        self.dataDense = A 
        self.L = np.linalg.cholesky(A) 
    def solve(self,x):
        z = np.linalg.solve(self.L,x)
        return np.linalg.solve(self.L.transpose(),z)        
###############################################################################       
class KPLUS:
    def __init__(self,Areg):
        self.__iAreg = spla.splu(Areg)  
#                
    def __mul__(self,b):
        return self.__iAreg.solve(b) 
###############################################################################       
class  KPLUS_HTFETI:
    def __init__(self,Kplus,B0,R,S_alpha): 
        self.Kplus  = Kplus
        self.B0     = B0
        self.R      = R
        for j in range(len(B0)):  
            if (j==0):       
                self.G0 = -sparse.csc_matrix.dot(B0[j],R[j])   
            else:
                G0i = -sparse.csc_matrix.dot(B0[j],R[j])    
                self.G0      = sparse.hstack((self.G0,G0i))  
        self.G0 = self.G0.tocsc()
        self.B0_Kplus = []
#                
        F0 = 0 
        for j in range(len(B0)): 
            B0_array = B0[j].toarray()
            self.B0_Kplus.append(np.zeros(B0_array.shape))
            for i in range(B0_array.shape[0]):
                self.B0_Kplus[j][i,:] = self.Kplus[j]*B0_array[i,:]
            F0 = F0 + np.dot(self.B0_Kplus[j],B0_array.transpose())
        self.iF0 = DENS_SOLVE(F0)
#
        S_alpha_from_ESPRESO = False
        if (not S_alpha_from_ESPRESO):
            G0d = self.G0.toarray();iF0_G0d = self.iF0.solve(G0d)        
            S_alpha = np.dot(G0d.T,iF0_G0d)
            rho = S_alpha[0,0]
            for i in range(R[0].shape[1]):
                S_alpha[i,i] += rho
#
        self.iS_alpha = DENS_SOLVE(S_alpha)           
#     
    def __mul__(self,f):
#
        d0 = np.zeros(self.B0[0].shape[0])        
        for i in range(len(self.Kplus)):
            if i == 0:
                e0 = sparse.csc_matrix.dot(self.R[i].transpose(),-f[i])
            else:
                e0i = sparse.csc_matrix.dot(self.R[i].transpose(),-f[i])
                e0 = np.concatenate((e0,e0i))
            d0 += np.dot(self.B0_Kplus[i],f[i])
#        
        tmpV = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))-e0
        alpha0 = self.iS_alpha.solve(tmpV)
        lam0 = self.iF0.solve(d0-sparse.csc_matrix.dot(self.G0,alpha0)) 
#
        cnt = 0
        uu = []
        for j in range(len(self.Kplus)):
            fj_B0t_lam0 = f[j]-sparse.csc_matrix.dot(self.B0[j].transpose(),lam0)
            uu.append(self.Kplus[j]*(fj_B0t_lam0))
            ind0 = np.arange(0,self.R[j].shape[1])+cnt            
            uu[j] += sparse.csc_matrix.dot(self.R[j],alpha0[ind0])
            cnt += self.R[j].shape[1]            
        return uu
###############################################################################      
class FETIOPERATOR: 
    def __init__(self,Kplus,B):
        self.Kplus  = Kplus
        self.B      = B        
    def __mul__(self,x_in):
        x_out = np.zeros(self.B[0][0].shape[0])
        for i in range(len(self.Kplus)):
            for j in range(len(self.Kplus[i])):
                x = sparse.csc_matrix.dot(self.B[i][j].transpose(),x_in)
                x = self.Kplus[i][j]*x
                x_out += sparse.csc_matrix.dot(self.B[i][j],x)
        return x_out
###############################################################################
class FETIOPERATOR_HTFETI:
    def __init__(self,Kplus_HTFETI,B1):
        self.Kplus_HTFETI = Kplus_HTFETI
        self.B1           = B1
        
    def __mul__(self,x_in):
        x_out = np.zeros(self.B1[0][0].shape[0])
        for i in range(len(self.Kplus_HTFETI)):
            xx = []
            for j in range(len(self.B1[i])):
                xx.append(sparse.csc_matrix.dot(self.B1[i][j].transpose(),x_in)) 
            x = self.Kplus_HTFETI[i]*xx
            for j in range(len(self.B1[i])):
                x_out += sparse.csc_matrix.dot(self.B1[i][j],x[j])           
#
        return x_out        
###############################################################################       
class COARSE_PROBLEM:
    def __init__(self,B,R):   
        for i in range(len(B)):
            for j in range(len(B[i])):  
                if (i==0 and j==0):
                    self.G = -sparse.csc_matrix.dot(B[i][j],R[i][j])   
                else:
                    Gi = -sparse.csc_matrix.dot(B[i][j],R[i][j])    
                    self.G      = sparse.hstack((self.G,Gi))  
        self.G = self.G.tocsc()
###############################################################################       
class COARSE_PROBLEM_HTFETI:
    def __init__(self,B1,R):   
        for i in range(len(B1)):
            for j in range(len(B1[i])):
                if (j==0):
                    Gj = sparse.csc_matrix.dot(B1[i][j],-R[i][j])   
                else:
                    Gj = Gj + sparse.csc_matrix.dot(B1[i][j],-R[i][j])    
            if i==0:                
                self.G      = Gj
            else:
                self.G      = sparse.hstack((self.G,Gj))                
        self.G = self.G.tocsc()
###############################################################################       
class PROJ:
    def __init__(self,G,iGtG):
        self.G = G
        self.iGtG = iGtG
    def __mul__(self,x):
        y = x - sparse.csc_matrix.dot(self.G,\
                (self.iGtG.solve(sparse.csc_matrix.dot(self.G.transpose(),x))))
        return y        
###############################################################################
class PREC_DIR_OR_LUMPED:
    def __init__(self,K,Schur,B,weight,index_weight):
        self.K  = K
        self.Schur  = Schur
        self.B  = B
        self.w  = weight
        self.iw = index_weight
#        
    def __mul__(self,x_in):
        if conf.precondDualSystem!='none':    
            x_out = np.zeros(self.B[0][0].shape[0])
            for i in range(len(self.K)):
                for j in range(len(self.K[i])):
                    x0  = x_in.copy()
                    x0[self.iw[i][j]] *= self.w[i][j]
                    x   = sparse.csc_matrix.dot(self.B[i][j].transpose(),x0)
                    if conf.precondDualSystem=='dirichlet':
                        x[self.Schur[i][j][0]]  =  \
                        np.dot(self.Schur[i][j][1],x[self.Schur[i][j][0]])
                    else:
                        x   = sparse.csc_matrix.dot(self.K[i][j],x)
                        
                    x1  = sparse.csc_matrix.dot(self.B[i][j],x)
                    x_out[self.iw[i][j]]+=x1[self.iw[i][j]]*self.w[i][j]
        else:
            x_out = x_in.copy() 
        return x_out        
###############################################################################
def pcgp(F, d, G, e, Prec, eps0, maxIt,disp,graph):
#
#
    GtG     = sparse.csc_matrix.dot(G.transpose(),G).tocsc()
    iGtG    = spla.splu(GtG)
    Proj    = PROJ(G,iGtG)
    lamIm   = sparse.csc_matrix.dot(G,iGtG.solve(e))
    nDual   = lamIm.shape[0]
    lam     = lamIm.copy()
    g       = F*lam - d     
    Pg      = Proj*g 
    MPg     = Prec*Pg
    PMPg    = Proj*MPg
#
    sqrt_gtPMPg0 = np.sqrt(np.dot(g,PMPg))
#
    if (np.dot(g,PMPg)<0): 
        raise SystemExit("Problem, precond. M is unsymmetric. Change it.") 
    sqrt_gtPg0      = np.sqrt(np.dot(g,Pg)) 
    vec_normed_g    = np.zeros(nDual) 
    vec_staggnat    = np.zeros(nDual)
#
    w     = PMPg.copy()
    strFormat= '%3.5f'
#
    if disp: 
        print('sqrt_gtPg0: %3.5e' %   (sqrt_gtPg0))
        print('sqrt_gtPMPg0:  %3.5e' % sqrt_gtPMPg0)
        print('    i      |r|        r         e         stagnation ') 
        print('%5d   %3.5f   %3.3e   %3.6f     %3.5f' % (1,1,sqrt_gtPg0,eps0,-1))     
#   
    for i in range(nDual):
#        
        Fw          = F*w
        rho         = -np.dot(g,PMPg)/np.dot(w,Fw)
        lam         += w * rho        
        gprev       = g.copy()                  
#         
        PMPgprev    = PMPg.copy()      
        g           += Fw * rho
        Pg           = Proj*g
        PMPg        = Proj*(Prec*Pg)

        gtPMPg      = np.dot(g,PMPg)
        gamma       = gtPMPg/np.dot(gprev,PMPgprev)
        w           = PMPg + w * gamma 
# 
        if (np.dot(g,PMPg)<0): 
            raise SystemExit("Problem, precond. M is unsymmetric. Change it.") 
#  
        sqrt_gtPg = np.sqrt(np.dot(g,Pg))
        normed_gi   = sqrt_gtPg/sqrt_gtPg0

        vec_normed_g[i]    = normed_gi
#       
        is_stagnating = np.log2(normed_gi)/(i+2)
        vec_staggnat[i]    = is_stagnating
#
        if np.log10(normed_gi/vec_normed_g[:i+1].min()) > 2:
            print('... stagnate',end='')
            break
#
        if disp:
            print('%5d   %3.5f   %3.3e   %3.6f     %3.5f' % \
                            (i+2,normed_gi,sqrt_gtPg,eps0,is_stagnating))             
#            
        if normed_gi<eps0:
            break
        if i==maxIt:
            print('PCPG does not converge within maxNumbIter (',
                       conf.maxIt_dual_feti,').')
            break
    alpha = iGtG.solve(sparse.csc_matrix.dot(G.transpose(),d-F*lam))
    numbOfIter = i
#   
    if graph:
        plt.subplot(2,1,1)    
        plt.plot(vec_staggnat[:i]);plt.ylim([vec_staggnat[:i].min()*2,0])
        plt.subplot(2,1,2)
        plt.semilogy(vec_normed_g[:i])#;plt.ylim([vec_staggnat[:i].min()*2,0])   
#  
    return lam, alpha, numbOfIter
###############################################################################       
def feti(K,Kreg,f,Schur,B,weight,index_weight,R):
#        
    maxIt   = conf.maxIt_dual_feti
    eps0    = conf.eps_dual_feti   
#              
    CP      = COARSE_PROBLEM(B,R)   
    d       = np.zeros(B[0][0].shape[0])    
    Kplus   = []
#    
    for i in range(len(K)):
        Kplus.append([])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS(Kreg[i][j]))
            if (i==0 and j==0):
                e = sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])
            else:
                e = np.concatenate((e,sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])))
            d += sparse.csc_matrix.dot(B[i][j],Kplus[i][j]*f[i][j])
#
    F       = FETIOPERATOR(Kplus,B)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B,weight,index_weight)
#     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
#
    uu = []
    cnt = 0
    print('size(lam):',lam.shape)
    print('size(B):',B[0][0].shape)
#
    delta = 0.0
    norm_f = 0.0
    for i in range(len(K)):
        uu.append([])
        for j in range(len(K[i])):
            ind = np.arange(0,R[i][j].shape[1]) + cnt
            f_BtLam_i_j = f[i][j]-sparse.csc_matrix.dot(B[i][j].transpose(),lam)
            KplusBtLam_i_j=Kplus[i][j]*f_BtLam_i_j
            R_alpha_i_j = sparse.csc_matrix.dot(R[i][j],alpha[ind])
            uu[i].append(KplusBtLam_i_j+R_alpha_i_j)
            cnt += R[i][j].shape[1]
            delta += np.linalg.norm(sparse.csc_matrix.dot(K[i][j],uu[i][j])-f_BtLam_i_j)                
            norm_f += np.linalg.norm(f[i][j])
#     
    print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam
###############################################################################    
def hfeti(K,Kreg,f,Schur,B0,B1,weight,index_weight,R,mat_S0):
#        
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
#                           
    CP  = COARSE_PROBLEM_HTFETI(B1,R)   
    d   = np.zeros(B1[0][0].shape[0])    
#   
    Kplus = []
    for i in range(len(K)):
        Kplus.append([])
        e_i_j = np.zeros(R[i][0].shape[1])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS(Kreg[i][j]))           
            e_i_j += sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])
        if (i==0):
            e = e_i_j
        else:
            e = np.concatenate((e,e_i_j))
#           
    Kplus_HTFETI = []       
    for i in range(len(K)):
        Kplus_HTFETI.append(KPLUS_HTFETI(Kplus[i],B0[i],R[i],mat_S0[i]))
#
    for i in range(len(K)):
        Kpl_ = Kplus_HTFETI[i]*f[i]
        for j in range(len(K[i])):
            d += sparse.csc_matrix.dot(B1[i][j],Kpl_[j])
#    
    F       = FETIOPERATOR_HTFETI(Kplus_HTFETI,B1)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B1,weight,index_weight)
     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
#    
    uu = []
    cnt = 0
    print('size(lam):',lam.shape)
    print('size(B):',B1[0][0].shape)
    print('type(B):',type(B1[0][0]))
    delta = 0.0
    norm_f = 0.0
    Kplus_f_B1t_lam = []
    cnt = 0
    for i in range(len(K)):
        B1t_lam = []
        for j in range(len(K[i])):
            B1t_lam.append(f[i][j]-sparse.csc_matrix.dot(B1[i][j].transpose(),lam))
        Kplus_f_B1t_lam.append(Kplus_HTFETI[i]*B1t_lam)
    uu = []   
    for i in range(len(K)):
        uu.append([])
        ind = np.arange(0,R[i][0].shape[1]) + cnt
        for j in range(len(K[i])):    
            R_alpha = sparse.csc_matrix.dot(R[i][j],alpha[ind])
            uu[i].append(Kplus_f_B1t_lam[i][j]+R_alpha)
        cnt += R[i][0].shape[1]
#    for i in range(len(K)):
#        uu.append([])
#        for j in range(len(K[i])):
#            ind = np.arange(0,R[i][j].shape[1]) + cnt
#            f_BtLam_i_j = f[i][j]-sparse.csc_matrix.dot(B1[i][j].transpose(),lam)
#            KplusBtLam_i_j=Kplus[i][j]*f_BtLam_i_j
#            R_alpha_i_j = sparse.csc_matrix.dot(R[i][j],alpha[ind])
#            uu[i].append(KplusBtLam_i_j+R_alpha_i_j)
#            cnt += R[i][j].shape[1]
#            delta += np.linalg.norm(sparse.csc_matrix.dot(K[i][j],uu[i][j])-f_BtLam_i_j)                
#            norm_f += np.linalg.norm(f[i][j])
        
#    u = Kplus_sub*f_m_BtLam + Roperator.mult(alpha)      
#     
#    delta = np.linalg.norm(Koperator.mult(u)-f_m_BtLam)
#    normf = np.linalg.norm(f)
#    print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam
############################################################################### 