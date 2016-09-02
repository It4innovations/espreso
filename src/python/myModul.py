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

class DENS_SOLVE_FOR_LU:   
    # solving of A*x = b with A = L*Lt
    # x = Lt^-1 * ( L^-1 * b)
    def __init__(self, A): 
        self.iA = spla.splu(sparse.csc_matrix( A )) 
    def solve(self,x): 
        return self.iA.solve(x)   
   
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
    def __init__(self,K,Kplus,Kplust,B0,R1,R2,S_alpha): 
        self.K      = K
        self.Kplus  = Kplus
        self.Kplust = Kplust
        self.B0     = B0
        self.R1     = R1
        self.R2     = R2
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',  R1[0].shape)

        for j in range(len(B0)):  
            if (j==0):       
                self.G01 = -sparse.csc_matrix.dot(B0[j],R1[j])  
                self.G02 = -sparse.csc_matrix.dot(B0[j],R2[j])   
            else:
                G01i = -sparse.csc_matrix.dot(B0[j],R1[j])  
                G02i = -sparse.csc_matrix.dot(B0[j],R2[j])    
                self.G01      = sparse.hstack((self.G01,G01i))   
                self.G02      = sparse.hstack((self.G02,G02i))  
        self.G01 = self.G01.tocsc()
        self.G02 = self.G02.tocsc()
        self.B0_Kplus = []
        self.B0_Kplust = []
#                
       
        
        
        iF0 = np.array([0],dtype=np.int32)
        jF0 = np.array([0],dtype=np.int32)
        vF0 = np.array([0])
        for j in range(len(B0)):  
            if False:
                B0_array = B0[j].toarray()
                self.B0_Kplus.append(np.zeros(B0_array.shape))
                self.B0_Kplust.append(np.zeros(B0_array.shape))
                for i in range(B0_array.shape[0]):
                    if B0_array[i,:].any():
                        self.B0_Kplus[j][i,:] = self.Kplust[j]*B0_array[i,:]  
                        self.B0_Kplust[j][i,:] = self.Kplus[j]*B0_array[i,:]  
                F0i = np.dot(self.B0_Kplus[j],B0_array.transpose())  
                F0 += F0i
            else:
                B0_bool = B0[j]!=0
                _lind = B0_bool.sum(1) 
                _ind = _lind.nonzero()[0].A1
                _n = _ind.shape[0]
                self.B0_Kplus.append(np.zeros(B0[j].shape))
                self.B0_Kplust.append(np.zeros(B0[j].shape))
                _j = np.zeros(_n*_n)
                for i in range(_ind.shape[0]):
                    B0_i = B0[j].getrow(_ind[i]).toarray()[0]  
                    self.B0_Kplus[j][_ind[i],:] = self.Kplust[j]*B0_i  
                    self.B0_Kplust[j][_ind[i],:] = self.Kplus[j]*B0_i  
                    _j[i*_n:(i+1)*_n]=_ind
                     
                _B0 = B0[j].tocsr()[_ind,:].tocsc() 
                F0i = sparse.csc_matrix.dot(_B0,self.B0_Kplust[j][_ind,:].T) 
                
                iF0 = np.concatenate((iF0,_ind.repeat(_n)))  
                jF0 = np.concatenate((jF0,_j))               
                vF0 = np.concatenate((vF0,F0i.ravel()))                
              
        _m = B0[0].shape[0]
        F0_sp = sparse.csc_matrix((vF0,(iF0,jF0)),shape=(_m,_m)).tocsc()
#            
        self.iF0  = spla.splu(F0_sp)
        self.iF0t = spla.splu(F0_sp.T)
#                
#        np.savetxt('F0_sp',F0_sp.toarray())        
        
        S_alpha_from_ESPRESO = False
#         
        if (not S_alpha_from_ESPRESO):
#           
            G01d = self.G01.toarray()
            G02d = self.G02.toarray()
            iF0_G01d = np.zeros(G01d.shape)
#            
            for i in range(G01d.shape[1]):
                iF0_G01d[:,i] = self.iF0.solve(G01d[:,i])        
            S_alpha = np.dot(G02d.T,iF0_G01d)
#            
            US, sS, VS = np.linalg.svd(S_alpha)
#            
            Hl = -US[:,-1]
            Hr = -VS[-1,:]
# #  ##############         
            LAMr = self.iF0.solve(sparse.csc_matrix.dot(self.G01,Hr))
            
            
            ttt = sparse.csc_matrix.dot(self.G02,Hl)            
            LAMl = self.iF0t.solve(ttt)
#            
            
            #np.savetxt('LAMr',LAMr) 
            #np.savetxt('LAMl',LAMl) 
            self.R1b = list(R1)
            self.R2b = list(R2) 
#            
            for j in range(len(self.R1b)):           
#                
                tmp1 = \
                    (self.Kplus[j]*(sparse.csc_matrix.dot(B0[j].T,LAMr)) + \
                    R1[j].data*Hr[j])
                ii = np.arange(0,tmp1.shape[0])
                jj = np.zeros(tmp1.shape[0]) 
#                
                self.R1b[j] = sparse.csc_matrix((tmp1,(ii,jj)),shape=(tmp1.shape[0],1)).tocsc()
#
                
#                B0tLAMl = sparse.csc_matrix.dot(B0[j].T,LAMl)
#                Kplust_B0tLAMl = self.Kplust[j] * B0tLAMl 
#                R2_Hl = R2[j].toarray()*Hl[j]
#                tmpp = -Kplust_B0tLAMl + R2_Hl
                vvv = sparse.csc_matrix.dot(B0[j].T,LAMl)
                Kplust_vvv = self.Kplust[j]*(vvv)
                _R_H = (R2[j].data*Hl[j])                
                tmp2 =  (Kplust_vvv +  _R_H)
#                    
                self.R2b[j] = sparse.csc_matrix((tmp2,(ii,jj)),shape=(tmp2.shape[0],1)).tocsc()
#
            print('LAMl=',LAMl)
            for j in range(len(self.R1b)):
#                
                vvr = self.K[j]*self.R1b[j].data + \
                        sparse.csc_matrix.dot(B0[j].T,LAMr)
#                
                vvl = self.K[j].T*self.R2b[j].data + \
                        sparse.csc_matrix.dot(B0[j].T,LAMl)
#
                print('a||K*Rr  + B0t*Lamr|| = ',np.linalg.norm(vvr))
                print('b||Kt*Rl + B0t*Laml|| = ',np.linalg.norm(vvl))
#            
            rho = S_alpha[0,0]
            for i in range(R1[0].shape[1]):
                S_alpha[i,i] += rho
#
# 
        np.savetxt('S_alpha',S_alpha)
        self.iS_alpha = DENS_SOLVE_FOR_LU(S_alpha)           
#     
    def __mul__(self,f):
#
        d0 = np.zeros(self.B0[0].shape[0])        
        for i in range(len(self.Kplus)):
            if i == 0:
                e0 = sparse.csc_matrix.dot(self.R2[i].transpose(),-f[i])
            else:
                e0i = sparse.csc_matrix.dot(self.R2[i].transpose(),-f[i])
                e0 = np.concatenate((e0,e0i))
            d0 += np.dot(self.B0_Kplus[i],f[i])
#        
        tmpV = sparse.csc_matrix.dot(self.G02.transpose(),self.iF0.solve(d0))-e0
        alpha0 = self.iS_alpha.solve(tmpV)
        lam0 = self.iF0.solve(d0-sparse.csc_matrix.dot(self.G01,alpha0)) 
#
        cnt = 0
        uu = []
        for j in range(len(self.Kplus)):
            fj_B0t_lam0 = f[j]-sparse.csc_matrix.dot(self.B0[j].transpose(),lam0)
            uu.append(self.Kplus[j]*(fj_B0t_lam0))
            ind0 = np.arange(0,self.R1[j].shape[1])+cnt            
            uu[j] += sparse.csc_matrix.dot(self.R1[j].data,alpha0[ind0])
            cnt += self.R1[j].shape[1]            
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
    def __init__(self,G,iGtG):
        self.G = G
        self.iGtG = iGtG
    def __mul__(self,x):
        y = x - sparse.csc_matrix.dot(self.G,\
                (self.iGtG.solve(sparse.csc_matrix.dot(self.G.transpose(),x))))
        return y          
###############################################################################
class PF_OPERATOR:
    def __init__(self,F,Proj):
        self.F = F
        self.Proj = Proj

    def __mul__(self,x_in):
        return self.Proj*self.F*x_in     
        
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
    
def GMRes(A, b, x0, e, nmax_iter, restart=None):
    
    
    r = b -  A*x0 

    x = []
    q = [0] * (nmax_iter)

    #x.append(r)
    #x  = np.zeros(b.shape[0])
    q[0] = r / np.linalg.norm(r)

    h = np.zeros((nmax_iter + 1, nmax_iter))
    
    for k in range(nmax_iter):
        y = np.asarray(A*q[k]).reshape(-1)

        for j in range(k):
            h[j, k] = np.dot(q[j], y)
            y = y - h[j, k] * q[j]
        h[k + 1, k] = np.linalg.norm(y)
        if (h[k + 1, k] != 0 and k != nmax_iter - 1):
            q[k + 1] = y / h[k + 1, k]

        b = np.zeros(nmax_iter + 1)
        b[0] = np.linalg.norm(r)

        result = np.linalg.lstsq(h, b)[0]

        x = (np.dot(np.asarray(q).transpose(), result) + x0)
    niter = k
    return x, niter    
    
    
    
    
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
        return Proj*(F*(Proj*v))
    
    PF         = LinearOperator((nDual, nDual), matvec=mv,dtype=np.float64)
    PF2        = PF_OPERATOR(F,Proj)
    
    d_          = Proj * (d - F * lamIm)
    counter     = gmres_counter()
    if 1:    
        lamKer,out  = sparse.linalg.gmres(PF,d_,d_*0,1e-11) 
        niter = 1
    else:
        lamKer, niter = GMRes(PF, d_, d_*0, 1e-4, 433, restart=None)
    lam         = lamKer + lamIm
    print('it: ',niter)
    print('||lam||=',np.linalg.norm( lam) )
    
    print('lam: ')
    print(lam)
    
    
    if os.path.isfile('lam_old.txt'):     
        lam_old = np.loadtxt('lam_old.txt')
        if lam.shape[0] == lam_old.shape[0]:
            print('norm = ',np.linalg.norm(lam-lam_old)/np.linalg.norm(lam))
    np.savetxt('lam_old.txt',lam)
    
    ones = np.ones(lam.shape[0])
    Fones = F * ones
    F_lam = F*lam
    dmFlam = d - F_lam
    alpha = iG2tG1.solve(sparse.csc_matrix.dot(G2.transpose(),dmFlam))
    
    #print('-----------------------', iG2tG1.solve(np.array([1])))
    print('aaa-', sparse.csc_matrix.dot(G2.transpose(),dmFlam))
    
    print('alpha = ', alpha)
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
    d        = np.zeros(B[0][0].shape[0])    
    dc       = np.zeros(B[0][0].shape[0])    
    Kplus    = []
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
def hfeti_unsym(K,Kreg,f,Schur,B0,B1,c,weight,index_weight,R1,R2,mat_S0):
#        
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
# 
    d           = np.zeros(B1[0][0].shape[0])    
    dc          = np.zeros(B1[0][0].shape[0])    
#    
    Kplus  = []
    Kplust = []
    for i in range(len(K)):
        Kplus.append([])
        Kplust.append([])
# 
        e_i_j = np.zeros(R2[i][0].shape[1])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS(Kreg[i][j]))
            tmpKt = Kreg[i][j].copy()
            Kplust[i].append(KPLUS(tmpKt.T))           
            e_i_j += sparse.csc_matrix.dot(R2[i][j].transpose(),-f[i][j])
        if (i==0):
            e = e_i_j
        else:
            e = np.concatenate((e,e_i_j))
#            
    for i in range(len(K)):
        for j in range(len(K[i])):
            print('norm(KR1):  ',np.linalg.norm(\
             sparse.csc_matrix.dot(K[i][j],R1[i][j].toarray())))
            print('norm(KtR2): ',np.linalg.norm(\
             sparse.csc_matrix.dot(K[i][j].T,R2[i][j].toarray())))
#             
    Kplus_HTFETI = []       
    for i in range(len(K)): 
        Kplus_HTFETI.append(KPLUS_HTFETI(K[i],Kplus[i],Kplust[i],B0[i],R1[i],R2[i],mat_S0[i]))
#        
        
        
    R1b = []
    R2b = []
    for i in range(len(K)):
        R1b.append([])
        R2b.append([])
        for j in range(len(K[i])):
#            print('a',R1b[i][j])
#            print('b',Kplus_HTFETI[i].R1b[j])
#            print('a',R2b[i][j])
#            print('b',Kplus_HTFETI[i].R2b[j])
            R1b[i].append(Kplus_HTFETI[i].R1b[j].copy())
            R2b[i].append(Kplus_HTFETI[i].R2b[j].copy())
#            print('c',R2b[i][j])
            
            
    CP1         = COARSE_PROBLEM_HTFETI(B1,R1b)  
    CP2         = COARSE_PROBLEM_HTFETI(B1,R2b) 
#
    for i in range(len(K)):
        Kpl_ = Kplus_HTFETI[i]*f[i]
        for j in range(len(K[i])):
            d += sparse.csc_matrix.dot(B1[i][j],Kpl_[j])
            dc[index_weight[i][j]] = c[i][j]
#    
    print('c = ', dc)
    d      -= dc
    F       = FETIOPERATOR_HTFETI(Kplus_HTFETI,B1)
    Prec    = PREC_DIR_OR_LUMPED(K,Schur,B1,weight,index_weight)
     
    #lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
    #Prec    = PREC_DIR_OR_LUMPED(K,Schur,B1,weight,index_weight)
#     
    lam, alpha, numbOfIter = gmres(F,d, CP1.G,CP2.G, e, Prec,eps0,maxIt,True,False)

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
        ind = np.arange(0,R1[i][0].shape[1]) + cnt
        for j in range(len(K[i])):    
            R_alpha = sparse.csc_matrix.dot(R1b[i][j],alpha[ind])
            uu[i].append(Kplus_f_B1t_lam[i][j]+R_alpha)
        cnt += R1[i][0].shape[1]


    return uu,lam
###############################################################################  
def hfeti(K,Kreg,f,Schur,B0,B1,c,weight,index_weight,R1,R2,mat_S0):
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
        
        e_i_j = np.zeros(R2[i][0].shape[1])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS(Kreg[i][j]))           
            e_i_j += sparse.csc_matrix.dot(R2[i][j].transpose(),-f[i][j])
        if (i==0):
            e = e_i_j
        else:
            e = np.concatenate((e,e_i_j))
#           
    Kplus_HTFETI = []       
    for i in range(len(K)):
        Kplus_HTFETI.append(KPLUS_HTFETI(Kplus[i],B0[i],R1[i],R2[i],mat_S0[i]))
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
