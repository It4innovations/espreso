import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import config_espreso_python as conf
import pylab as plt
import os.path

#import multiprocessing

#import sys  
#
#
def load_matrix(path,str0,i,j,makeSparse,makeSymmetric): 
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    print(pathToFile)
    isFile = os.path.isfile(pathToFile)
    
    if isFile:
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
    else:
      tmp = []
#
    return tmp
###############################################################################
def load_vector(path,str0,i,j):
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    print(pathToFile)
    isFile = os.path.isfile(pathToFile)
    if isFile:
      tmp = np.loadtxt(pathToFile)
    else:
      tmp = []
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
class SPARSE_SOLVE: 
    # solving of A*x = b with A = L*Lt
    # x = Lt^-1 * ( L^-1 * b)
    def __init__(self, A):
        self.A = A 
        self.dA = np.sqrt(A.diagonal());
        nA = A.shape[0]
        S  = sparse.spdiags(1.0/self.dA,0,nA,nA)
        SAS = sparse.csc_matrix.dot(S,A)
        SAS = sparse.csc_matrix.dot(SAS,S)
        self.iA = spla.splu(SAS) 
    def solve(self,b):
        z = self.iA.solve(b)*self.dA
        return z        
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
        
        iF0 = np.array([0],dtype=np.int32)
        jF0 = np.array([0],dtype=np.int32)
        vF0 = np.array([0])
        for j in range(len(B0)):  
            if False:
                B0_array = B0[j].toarray()
                self.B0_Kplus.append(np.zeros(B0_array.shape))
                for i in range(B0_array.shape[0]):
                    if B0_array[i,:].any():
                        self.B0_Kplus[j][i,:] = self.Kplus[j]*B0_array[i,:]  
                F0i = np.dot(self.B0_Kplus[j],B0_array.transpose())  
                F0 += F0i
            else:
                B0_bool = B0[j]!=0
                _lind = B0_bool.sum(1) 
                _ind = _lind.nonzero()[0].A1
                _n = _ind.shape[0]
                self.B0_Kplus.append(np.zeros(B0[j].shape))
                _j = np.zeros(_n*_n)
                for i in range(_ind.shape[0]):
                    B0_i = B0[j].getrow(_ind[i]).toarray()[0]  
                    self.B0_Kplus[j][_ind[i],:] = self.Kplus[j]*B0_i  
                    _j[i*_n:(i+1)*_n]=_ind
                     
                _B0 = B0[j].tocsr()[_ind,:].tocsc() 
                F0i = sparse.csc_matrix.dot(_B0,self.B0_Kplus[j][_ind,:].T) 
                
                iF0 = np.concatenate((iF0,_ind.repeat(_n)))  
                jF0 = np.concatenate((jF0,_j))               
                vF0 = np.concatenate((vF0,F0i.ravel()))                
                
                
                #F0[np.ix_(_ind,_ind)] += F0i
              
        _m = B0[0].shape[0]
        F0_sp = sparse.csc_matrix((vF0,(iF0,jF0)),shape=(_m,_m)).tocsc()
            
        #np.savetxt('F0',F0)
#        np.savetxt('F0_sp',F0_sp.toarray())

#        self.iF0 = DENS_SOLVE(F0_sp.toarray())
#        self.iF0sp = SPARSE_SOLVE(F0_sp)
        self.iF0 = spla.splu(F0_sp)
#
#        ee  = np.ones(_m)
#        print(self.iF0.solve(ee))
#        print(self.iF0sp_.solve(ee))
#        err = np.linalg.norm(self.iF0sp_.solve(ee)-self.iF0.solve(ee))        
#        print(err)
        
        
        S_alpha_from_ESPRESO = False
        if (not S_alpha_from_ESPRESO):
            G0d = self.G0.toarray()
            iF0_G0d = np.zeros(G0d.shape)
            for i in range(G0d.shape[1]):
                iF0_G0d[:,i] = self.iF0.solve(G0d[:,i])        
            S_alpha = np.dot(G0d.T,iF0_G0d)
            rho = S_alpha[0,0]
            for i in range(R[0].shape[1]):
                S_alpha[i,i] += rho
#
#        np.savetxt('F0',F0)
        np.savetxt('S_alpha',S_alpha)
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
class PROJ_UNSYM:
    def __init__(self,G1,G2,iG2tG1):
        self.G1 = G1
        self.G2 = G2
        self.iG2tG1 = iG2tG1
    def __mul__(self,x):
        y = x - sparse.csc_matrix.dot(self.G1,\
                (self.iG2tG1.solve(sparse.csc_matrix.dot(self.G2.transpose(),x))))
        return y                  
###############################################################################
class PROJ:
    def __init__(self,G1,G2,iG1tG):
        self.G = G
        self.iGtG = iGtG
    def __mul__(self,x):
        y = x - sparse.csc_matrix.dot(self.G,\
                (self.iGtG.solve(sparse.csc_matrix.dot(self.G.transpose(),x))))
        return y          
###############################################################################
class PFP_OPERATOR:
    def __init__(self,F,Proj):
        self.F = F
        self.Proj = Proj

    def __mul__(self,x_in):
        return self.Proj*self.F*self.Proj*x_in     
        
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
class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
    def __call__(self, rk=None):
        self.niter += 1
        if self._disp:
            print('iter %3i\trk = %s' % (self.niter, str(rk)))    
    
    
    
    
    
    
def gmres(F, d, G1, G2, e, Prec, eps0, maxIt,disp,graph):
#
#
    from scipy.sparse.linalg import LinearOperator
    G2tG1       = sparse.csc_matrix.dot(G2.transpose(),G1).tocsc()
    iG2tG1      = spla.splu(G2tG1)
    Proj        = PROJ_UNSYM(G1,G2,iG2tG1)
    lamIm       = sparse.csc_matrix.dot(G1,iG2tG1.solve(e))
    nDual       = lamIm.shape[0]
    lam         = lamIm.copy()
    g           = F*lam - d     
    i           = 0
    #PFP         = PFP_OPERATOR(F, Proj)
    
    def mv(v):
        return Proj*(F*(v))
    
    PF         = LinearOperator((nDual, nDual), matvec=mv,dtype=np.float64)

    
    d_          = Proj * (d - F * lamIm)
    counter     = gmres_counter()
    lamKer,out  = sparse.linalg.gmres(PF,d_) 
    lam         = lamKer + lamIm
    print('out: ',out)
    print('it: ',counter.niter)
    print('||lam||=',np.linalg.norm( lam) )
    
    print('lam: ')
    print(lam-11)
    
    
    if os.path.isfile('lam_old.txt'):     
        lam_old = np.loadtxt('lam_old.txt')
        if lam.shape[0] == lam_old.shape[0]:
            print('norm = ',np.linalg.norm(lam-lam_old)/np.linalg.norm(lam))
    np.savetxt('lam_old.txt',lam)
    

    alpha = iG2tG1.solve(sparse.csc_matrix.dot(G2.transpose(),d-F*lam))
    numbOfIter = i
#   
#    if graph:
#        plt.subplot(2,1,1)    
#        plt.plot(vec_staggnat[:i]);plt.ylim([vec_staggnat[:i].min()*2,0])
#        plt.subplot(2,1,2)
#        plt.semilogy(vec_normed_g[:i])#;plt.ylim([vec_staggnat[:i].min()*2,0])   
#  
    return lam, alpha, numbOfIter
###############################################################################       
def feti_unsym(K,Kreg,f,Schur,B,c,weight,index_weight,R1,R2):
#        
    maxIt   = conf.maxIt_dual_feti
    eps0    = conf.eps_dual_feti   
#              
    CP1      = COARSE_PROBLEM(B,R1)  
    CP2      = COARSE_PROBLEM(B,R2)   
    d       = np.zeros(B[0][0].shape[0])    
    dc       = np.zeros(B[0][0].shape[0])    
    Kplus   = []
#    
    
    for i in range(len(K)):
        Kplus.append([])
        
        
        for j in range(len(K[i])):

            Kplus[i].append(KPLUS(Kreg[i][j]))
            if (i==0 and j==0):
                e = sparse.csc_matrix.dot(R1[i][j].transpose(),-f[i][j])
            else:
                e = np.concatenate((e,sparse.csc_matrix.dot(R1[i][j].transpose(),-f[i][j])))
            d += sparse.csc_matrix.dot(B[i][j],Kplus[i][j]*f[i][j])
            dc[index_weight[i][j]] = c[i][j]
#
    
    d      -= dc    
    F       = FETIOPERATOR(Kplus,B)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B,weight,index_weight)
#     
    lam, alpha, numbOfIter = gmres(F,d, CP1.G,CP2.G, e, Prec,eps0,maxIt,True,False)        
#
    uu = []
    cnt = 0
    #print('size(lam):',lam.shape)
    #print('size(B):',B[0][0].shape)
#
    delta = 0.0
    norm_f = 0.0
    for i in range(len(K)):
        uu.append([])
        for j in range(len(K[i])):
            ind = np.arange(0,R1[i][j].shape[1]) + cnt
            f_BtLam_i_j = f[i][j]-sparse.csc_matrix.dot(B[i][j].transpose(),lam)
            KplusBtLam_i_j=Kplus[i][j]*f_BtLam_i_j
            R_alpha_i_j = sparse.csc_matrix.dot(R1[i][j],alpha[ind])
            uu[i].append(KplusBtLam_i_j+R_alpha_i_j)
            cnt += R1[i][j].shape[1]
            delta += np.linalg.norm(sparse.csc_matrix.dot(K[i][j],uu[i][j])-f_BtLam_i_j)                
            norm_f += np.linalg.norm(f[i][j])
#     
    if np.abs(norm_f)>1e-10:
        print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam


def feti(K,Kreg,f,Schur,B,c,weight,index_weight,R):
#        
    maxIt   = conf.maxIt_dual_feti
    eps0    = conf.eps_dual_feti   
#              
    CP      = COARSE_PROBLEM(B,R)   
    d       = np.zeros(B[0][0].shape[0])    
    dc       = np.zeros(B[0][0].shape[0])    
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
            dc[index_weight[i][j]] = c[i][j]
#
    d      -= dc        
    F       = FETIOPERATOR(Kplus,B)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B,weight,index_weight)
#     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
#
    uu = []
    cnt = 0
    #print('size(lam):',lam.shape)
    #print('size(B):',B[0][0].shape)
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
    if np.abs(norm_f)>1e-10:
        print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam
###############################################################################    
def hfeti(K,Kreg,f,Schur,B0,B1,c,weight,index_weight,R,mat_S0):
#        
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
#                           
    CP  = COARSE_PROBLEM_HTFETI(B1,R)   
    d   = np.zeros(B1[0][0].shape[0])    
    dc  = np.zeros(B1[0][0].shape[0])    
#   
#    pool = multiprocessing.Pool()

    Kplus = []
    for i in range(len(K)):
        Kplus.append([])
#        Kplus.append(pool.map(KPLUS,Kreg[i]))
        
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
            dc[index_weight[i][j]] = c[i][j]
#    
    d      -= dc
    F       = FETIOPERATOR_HTFETI(Kplus_HTFETI,B1)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B1,weight,index_weight)
     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
#    
    uu = []
    cnt = 0
    #print('size(lam):',lam.shape)
    #print('size(B):',B1[0][0].shape)
    #print('type(B):',type(B1[0][0]))
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
#        ind = np.arange(0,R[i][j].shape[1]) + cnt
#        for j in range(len(K[i])):
#            f_BtLam_i_j = f[i][j]-sparse.csc_matrix.dot(B1[i][j].transpose(),lam)
#            KplusBtLam_i_j=Kplus[i][j]*f_BtLam_i_j
#            R_alpha_i_j = sparse.csc_matrix.dot(R[i][j],alpha[ind])
#            uu[i].append(KplusBtLam_i_j+R_alpha_i_j)
#            cnt += R[i][j].shape[1]
#            delta += np.linalg.norm(sparse.csc_matrix.dot(K[i][j],uu[i][j])-f_BtLam_i_j)                
#            norm_f += np.linalg.norm(f[i][j])
        
    #u = Kplus_sub*f_m_BtLam + Roperator.mult(alpha)      
     
    #delta = np.linalg.norm(Koperator.mult(u)-f_m_BtLam)
    #norm_f = np.linalg.norm(f)
#    if np.abs(norm_f)>1e-10:
#        print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam
############################################################################### 
