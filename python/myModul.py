import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import config_espreso_python as conf
import pylab as plt
#import sys  
 
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
class PREC_PRIMAL_SYSTEM:
    def __init__(self,list_A):        
        if conf.precondPrimalSystem=='diag':
            self.iAdiag = list_A[1] 
            self.flag=1
        elif conf.precondPrimalSystem == 'LU_SP' or \
                        conf.precondPrimalSystem == 'ILU':
            self.LU_precond=list_A[1]
            self.flag=1
        elif conf.precondPrimalSystem == 'none':
            self.flag=0
    def mult(self,x):
        if conf.precondPrimalSystem=='diag':
            return self.iAdiag*x
        elif conf.precondPrimalSystem == 'LU_SP':        
            return self.LU_precond.solve(x.astype(np.float32)) 
        elif conf.precondPrimalSystem == 'ILU': 
            if conf.single_precision:
                return self.LU_precond.solve(x.astype(np.float32))
            else:
                return self.LU_precond.solve(x)                
        elif conf.precondPrimalSystem == 'none':
            return x
###############################################################################           
class MULT_SPARSE_MATRIX:
    def __init__(self,list0):
        self.A = list0[0]
        self.n_list = len(list0)
        if self.n_list==2:
            self.R = list0[1]
            self.rho = self.A.diagonal().mean()
            
    def mult(self,x):
        if self.n_list==1:
            return sparse.csc_matrix.dot(self.A,x)            
        else:
            y = sparse.csc_matrix.dot(self.A,x)
            y += self.rho * sparse.csc_matrix.dot(self.R,\
                    sparse.csc_matrix.dot(self.R.transpose(),x))
            return y
    
###############################################################################       
class KPLUS_NEW:
    def __init__(self,Areg):
        self.iAreg = spla.splu(Areg)  
                
    def __mul__(self,b):
        return self.iAreg.solve(b) 
        
class KPLUS:     
    # 
    def __init__(self, A, Areg, R): 
        
        self.Areg               = Areg    
        self.A                  = A
        
        if conf.iterative_Kplus:  
            self.iAreg              = np.ones(self.A.shape[0])
            self.R                  = R            
            #
            if conf.precondPrimalSystem == 'diag':
                if conf.precondFrom_Areg_orA:
                    self.iAreg = self.Areg.diagonal()**-1
                else:      
                    if conf.mult_Areg_or_A_RRt:
                        self.iAreg = self.A.diagonal()**-1
                    else:
                        self.iAreg = (self.A.diagonal())**-1 
            elif conf.precondPrimalSystem == 'none':
                    self.iAreg = 0
            else:
                if conf.precondFrom_Areg_orA:
                    self.A_Prec_copy        = Areg.copy()  

                else: 
                    self.A_Prec_copy        = A.copy() 
#                    if conf.precondPrimalSystem != 'ILU':
#                        print('precondPrimalSystem has to be set to ILU')
#                        print(' ... will be pre-set')
#                        conf.precondPrimalSystem='ILU'
                #  
                if (conf.single_precision):
                    self.A_Prec_copy = self.A_Prec_copy.astype(np.float32)
                #
                if conf.precondPrimalSystem=='ILU':
                    self.iAreg = spla.spilu(self.A_Prec_copy,1e-4)
                elif conf.precondPrimalSystem=='LU_SP':
                    self.iAreg = spla.splu(self.A_Prec_copy)

                
        else:
            self.iAreg=[] 
            for i in range(len(self.A)):
                self.iAreg.append([]) 
                for j in range(len(self.A[i])):
                    self.iAreg[i].append(spla.splu(self.Areg[i][j]))  

    def __mul__(self,b,use_subset=-1): 
        if conf.iterative_Kplus:            
            maxIt   = conf.maxIt_dual_feti
            x0      = b*0 #self.iAreg.solve(b.astype(np.float32))
            Prec    = PREC_PRIMAL_SYSTEM([self.A,self.iAreg])
            
            
            if conf.methodToImproveSolByIterMethod[:2] == 'cg':

                if conf.mult_Areg_or_A_RRt:
                    mult_A  = MULT_SPARSE_MATRIX([self.Areg])
                else:
                    mult_A  = MULT_SPARSE_MATRIX([self.Areg,self.R])            
                if conf.methodToImproveSolByIterMethod == 'cg_x':
                    x, norm_g, numbOfIter = cg(mult_A,b,x0,self.R,\
                            conf.eps_iter_Kplus,maxIt,Prec) 
                elif conf.methodToImproveSolByIterMethod == 'cg_dx':
                    db = b - sparse.csc_matrix.dot(mult_A,x0)
                    dx, norm_g, numbOfIter = cg(self.A,db,x0*0,self.R,\
                            conf.eps_iter_Kplus,maxIt,Prec)
                    x = x0 + dx
            elif conf.methodToImproveSolByIterMethod == 'pcg_x':
  
                mult_A  = MULT_SPARSE_MATRIX([self.A])
                ee      = np.zeros(self.R.shape[1])

                x, lam, numbOfIter = \
                        pcgp(mult_A, b, self.R, ee, Prec, \
                        conf.eps_iter_Kplus, maxIt, False,False)
            else:
                print('chose from: cg_x, cg_dx or pcg_x')
                    
            print(numbOfIter,',',end='')   
            print('%d,' % (numbOfIter),end='')
            return x            
        else:
            x_out = np.zeros(b.shape[0])
            cnt=0
            if use_subset==-1:
                indOfSub = np.arange(len(self.iAreg))
            else:
                indOfSub = np.arange(use_subset)
            for i in indOfSub:
                for j in range(len(self.iAreg[i])):
                    nA = self.A[i][j].shape[0]
                    ind = np.arange(0,nA)+cnt
                    x_out[ind] = self.iAreg[i][j].solve(b[ind]) 
                    cnt += nA
#            else:
#                for i in range(use_subset[0]):
#                    for j in range(len(self.iAreg[i])):
#                        nA = self.A[i][j].shape[0]
#                        ind = np.arange(0,nA)+cnt
#                        x_out[ind] = self.iAreg[i][j].solve(b[ind]) 
#                        cnt += nA
            return x_out  
###############################################################################  
            
class MULT_BLOCK_DIAG:
    def __init__(self,A):
        self.A = A
        self.n = 0
        self.m = 0
        for i in range(len(self.A)):
            for j in range(len(self.A[i])):
                self.n += A[i][j].shape[0]
                self.m += A[i][j].shape[1]
        
    def mult(self,b):
        x_out = np.zeros(self.n)
        cnt_n=0
        cnt_m=0
        for i in range(len(self.A)):
            for j in range(len(self.A[i])):
                n = self.A[i][j].shape[0]
                m = self.A[i][j].shape[1]
                ind_n = np.arange(0,n)+cnt_n
                ind_m = np.arange(0,m)+cnt_m
                x_out[ind_n] = sparse.csc_matrix.dot(self.A[i][j],b[ind_m])                
                cnt_n += n                 
                cnt_m += m    
        return x_out  
        
    def mult_t(self,b):
        x_out = np.zeros(self.m)
        cnt_n=0
        cnt_m=0
        for i in range(len(self.A)):
            for j in range(len(self.A[i])):
                n = self.A[i][j].shape[0]
                m = self.A[i][j].shape[1]
                ind_n = np.arange(0,n)+cnt_n
                ind_m = np.arange(0,m)+cnt_m
                x_out[ind_m] = sparse.csc_matrix.dot(self.A[i][j].transpose(),b[ind_n])                
                cnt_n += n                 
                cnt_m += m    
        return x_out 
            
class MULT_BLOCK_RECTANGLE:
    def __init__(self,B):
        self.B = B
        self.n = self.B[0][0].shape[0]
        
        self.m = 0 
        for i in range(len(self.B)):
            for j in range(len(self.B[i])):
                self.m += self.B[i][j].shape[1]        
        
        
        
    def mult(self,b):
        x_out = np.zeros(self.n)
        cnt=0
        for i in range(len(self.B)):
            for j in range(len(self.B[i])):
                m_i = self.B[i][j].shape[1]
                ind = np.arange(0,m_i)+cnt
                x_out += sparse.csc_matrix.dot(self.B[i][j],b[ind])                
                cnt += m_i    
        return x_out  

    def mult_t(self,b):
        x_out = np.zeros(self.m)
        cnt=0
        for i in range(len(self.B)):
            for j in range(len(self.B[i])):
                nB = self.B[i][j].shape[1]
                ind = np.arange(0,nB)+cnt
                x_out[ind] += sparse.csc_matrix.dot(self.B[i][j].T,b)                
                cnt += nB    
        return x_out  


#            
#class KPLUS_HFETI:
#    def __init__(self, B0,G0,Kplus_sub,R,iF0,iS0):
#        self.B0 = B0
#        self.G0 = G0
#        self.nB0 = B0.shape[0]
#        self.Kplus_sub = Kplus_sub
#        self.R = R
#        self.iF0 = iF0
#        self.iS0 = iS0
#        
#        self.B_Kplus_sub = np.zeros(B0.shape)
#        for i in range(len(self.B0)):
#            for j in range(len(self.B[i])): 
#                B0_array = B0[i][j].toarray()
#        
#        
#        
#        for i in range(B0.shape[0]):
#            self.B_Kplus_sub[i,:] = self.Kplus_sub.solve(B0_array[i,:])
#        
#    def solve(self,bc):
#        for i in range(G0.shape[0]):
#            b = bc[:-self.nB0:]
#            c = bc[-self.nB0::]      
#            Kplus_b = self.Kplus_sub.solve(b)
#            d0 = sparse.csc_matrix.dot(self.B0,Kplus_b)-c 
#            e0 = -sparse.csc_matrix.dot(self.R.transpose(),b)   
#            G0tiF0d0 = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))
#            beta = self.iS0.solve(G0tiF0d0-e0) 
#            mu = self.iF0.solve(d0 - sparse.csc_matrix.dot(self.G0,beta))    
#            x = Kplus_b-np.dot(self.B_Kplus_sub.transpose(),mu) + \
#               sparse.csc_matrix.dot(self.R,beta)  
##       x = Kplus_b-self.Kplus_sub.solve(sparse.csc_matrix.dot(self.B0.transpose(),mu)) + \
##            sparse.csc_matrix.dot(self.R,beta)  
#     
#        return np.concatenate((x,mu))        
#        
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
                
        
        F0 = 0 
        for j in range(len(B0)): 
            B0_array = B0[j].toarray()
            self.B0_Kplus.append(np.zeros(B0_array.shape))
            for i in range(B0_array.shape[0]):
                self.B0_Kplus[j][i,:] = self.Kplus[j]*B0_array[i,:]
            F0 = F0 + np.dot(self.B0_Kplus[j],B0_array.transpose())
        self.iF0 = DENS_SOLVE(F0)


        S_alpha_from_ESPRESO = False
        if (not S_alpha_from_ESPRESO):
            G0d = self.G0.toarray();iF0_G0d = self.iF0.solve(G0d)        
            S_alpha = np.dot(G0d.T,iF0_G0d)
            rho = S_alpha[0,0]
            for i in range(R[0].shape[1]):
                S_alpha[i,i] += rho
                
#        for i in range(B0):
#            np.savetxt('B0_Kplus'+int(i)+'.txt',self.B0_Kplus[i])
#        np.savetxt('F0.txt',F0)
        self.iS_alpha = DENS_SOLVE(S_alpha)   
        
     
    def __mul__(self,f):

        d0 = np.zeros(self.B0[0].shape[0])        
        for i in range(len(self.Kplus)):
            if i == 0:
                e0 = sparse.csc_matrix.dot(self.R[i].transpose(),-f[i])
            else:
                e0i = sparse.csc_matrix.dot(self.R[i].transpose(),-f[i])
                e0 = np.concatenate((e0,e0i))
            d0 += np.dot(self.B0_Kplus[i],f[i])
        
        
        tmpV = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))-e0
        alpha0 = self.iS_alpha.solve(tmpV)
        lam0 = self.iF0.solve(d0-sparse.csc_matrix.dot(self.G0,alpha0)) 

#        np.savetxt('alpha0.txt',alpha0)
#        raise


        cnt = 0
        uu = []
        for j in range(len(self.Kplus)):
            fj_B0t_lam0 = f[j]-sparse.csc_matrix.dot(self.B0[j].transpose(),lam0)
            uu.append(self.Kplus[j]*(fj_B0t_lam0))
            ind0 = np.arange(0,self.R[j].shape[1])+cnt            
            uu[j] += sparse.csc_matrix.dot(self.R[j],alpha0[ind0])
            cnt += self.R[j].shape[1]        
            
        return uu
        

#
#        def __mul__(self,b):
#            for j in range(self.G0.shape[0]):
#                Kplus_b = Kplus[j]*b
#                d0 = sparse.csc_matrix.dot(self.B0[i],Kplus_b) 
#                e0 = -sparse.csc_matrix.dot(self.R[i].transpose(),b)   
#                G0tiF0d0 = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))
#                beta = self.iS0.solve(G0tiF0d0-e0) 
#                mu = self.iF0.solve(d0 - sparse.csc_matrix.dot(self.G0,beta))    
#                x = Kplus_b-np.dot(self.B_Kplus_sub.transpose(),mu) + \
#                   sparse.csc_matrix.dot(self.R,beta)
#            return x_out
      
      
      
#    def __mult__(self,x_in):
#        for i in range(G0.shape[0]):
#            b = bc[:-self.nB0:]
#            c = bc[-self.nB0::]      
#            Kplus_b = self.Kplus_sub.solve(b)
#            d0 = sparse.csc_matrix.dot(self.B0,Kplus_b)-c 
#            e0 = -sparse.csc_matrix.dot(self.R.transpose(),b)   
#            G0tiF0d0 = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))
#            beta = self.iS0.solve(G0tiF0d0-e0) 
#            mu = self.iF0.solve(d0 - sparse.csc_matrix.dot(self.G0,beta))    
#            x = Kplus_b-np.dot(self.B_Kplus_sub.transpose(),mu) + \
#               sparse.csc_matrix.dot(self.R,beta)
#        return x_out
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
#            x = self.Kplus_HTFETI[i].mult(xx)             
            x = self.Kplus_HTFETI[i]*xx
            for j in range(len(self.B1[i])):
                x_out += sparse.csc_matrix.dot(self.B1[i][j],x[j])
                
                
            #Kpl_ = self.Kplus_HTFETI[i].mult(f[i])
            #for j in range(len(K[i])):
            #    d += sparse.csc_matrix.dot(B1[i][j],Kpl_[j])    
        
        return x_out        
        
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
    
 

   
#class COARSE_PROBLEM_HTFETI_CLUSTER:
#    def __init__(self,B0,R):
#        self.G0 = []
#        for i in range(len(B0)):
#            for j in range(len(B0[i])):  
#                if (j==0):
#                    self.G0.append(-sparse.csc_matrix.dot(B0[i][j],R[i][j]))   
#                else:
#                    G0i = -sparse.csc_matrix.dot(B0[i][j],R[i][j])    
#                    self.G0[i]      = sparse.hstack((self.G0[i],G0i))   

 

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
    def __init__(self,K,B):
        self.K = K
        self.B = B
    def __mul__(self,x):
        if True:    
            y = x.copy()   
        else:
            y = sparse.csc_matrix.dot(self.B.transpose(),x)
            y = sparse.csc_matrix.dot(self.K,y)
            y = sparse.csc_matrix.dot(self.B,y)
        return y 
###############################################################################
def load_matrix(path,str0,i,j,makeSparse,makeSymmetric): 
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile, ndmin=2)
    if (tmp.shape[0]==1):
        tmp = []
    else:
        n = np.int32(tmp[0,0])   
        m = np.int32(tmp[0,1])
        I = tmp[1::,0]-1;    J = tmp[1::,1]-1;    V = tmp[1::,2]
    
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
        
    return tmp
###############################################################################
def load_vector(path,str0,i,j):
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile)
    return tmp
###############################################################################  
def pcgp(F, d, G, e, Prec, eps0, maxIt,disp,graph):
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
        print('%5d   %3.5f   %3.3e   %3.6f     %3.5f' % (1,1,sqrt_gtPMPg0,eps0,-1))     
#   
    for i in range(nDual):
        
        Fw          = F*w
        rho         = -np.dot(g,PMPg)/np.dot(w,Fw)
        lam         += w * rho        
        gprev       = g.copy()                  
#
         
        PMPgprev    = PMPg.copy()      
        g           += Fw * rho
        PMPg        = Proj*(Prec*(Proj*g))

        gtPMPg      = np.dot(g,PMPg)
        gamma       = gtPMPg/np.dot(gprev,PMPgprev)
        w           = PMPg + w * gamma 
# 
        if (np.dot(g,PMPg)<0): 
            raise SystemExit("Problem, precond. M is unsymmetric. Change it.") 
#  
        sqrt_gtPMPg = np.sqrt(gtPMPg)
        normed_gi   = sqrt_gtPMPg/sqrt_gtPg0

        vec_normed_g[i]    = normed_gi
        
        is_stagnating = np.log2(normed_gi)/(i+2)
        vec_staggnat[i]    = is_stagnating
#
        if np.log10(normed_gi/vec_normed_g[:i+1].min()) > 2:
            print('... stagnate',end='')
            break

        if disp:
            print('%5d   %3.5f   %3.3e   %3.6f     %3.5f' % \
                            (i+2,normed_gi,sqrt_gtPMPg,eps0,is_stagnating))             
            
            
        if normed_gi<eps0:
            break
        if i==maxIt:
            print('PCPG does not converge within maxNumbIter (',
                       conf.maxIt_dual_feti,').')
            break
    alpha = iGtG.solve(sparse.csc_matrix.dot(G.transpose(),d-F*lam))
    numbOfIter = i
    
    if graph:
        plt.subplot(2,1,1)    
        plt.plot(vec_staggnat[:i]);plt.ylim([vec_staggnat[:i].min()*2,0])
        plt.subplot(2,1,2)
        plt.semilogy(vec_normed_g[:i])#;plt.ylim([vec_staggnat[:i].min()*2,0])
#   
#  
    return lam, alpha, numbOfIter
    
    
    
  
def feti(K,Kreg,f,B,R,weight):
        
    weight = 1
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
              
    CP      = COARSE_PROBLEM(B,R)   
    d = np.zeros(B[0][0].shape[0])    
    Kplus = []
    for i in range(len(K)):
        Kplus.append([])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS_NEW(Kreg[i][j]))
            if (i==0 and j==0):
                e = sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])
            else:
                e = np.concatenate((e,sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])))
            d += sparse.csc_matrix.dot(B[i][j],Kplus[i][j]*f[i][j])


    F       = FETIOPERATOR(Kplus,B)
    Prec    = PREC_DIR_OR_LUMPED(K,B)
     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,True)        
#    
##    for i in range(len(f))
    uu = []
    cnt = 0
    print('size(lam):',lam.shape)
    print('size(B):',B[0][0].shape)
    print('type(B):',type(B[0][0]))
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
def hfeti(K,Kreg,f,B0,B1,R,mat_S0,weight):
#def feti(K,Kreg,f,B,R,weight):
        
    weight = 1
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
                            
    CP  = COARSE_PROBLEM_HTFETI(B1,R)   
    d   = np.zeros(B1[0][0].shape[0])    
    
    
    
    
    
    Kplus = []
    for i in range(len(K)):
        Kplus.append([])
        e_i_j = np.zeros(R[i][0].shape[1])
        for j in range(len(K[i])):
            Kplus[i].append(KPLUS_NEW(Kreg[i][j]))           
            e_i_j += sparse.csc_matrix.dot(R[i][j].transpose(),-f[i][j])
        if (i==0):
            e = e_i_j
        else:
            e = np.concatenate((e,e_i_j))
            
    Kplus_HTFETI = []       
    for i in range(len(K)):
        Kplus_HTFETI.append(KPLUS_HTFETI(Kplus[i],B0[i],R[i],mat_S0[i]))


#    np.savetxt('G0.txt',Kplus_HTFETI[0].G0.todense())


    for i in range(len(K)):
        Kpl_ = Kplus_HTFETI[i]*f[i]
        for j in range(len(K[i])):
            d += sparse.csc_matrix.dot(B1[i][j],Kpl_[j])
     
            
             
 
    
    
    F       = FETIOPERATOR_HTFETI(Kplus_HTFETI,B1)
    Prec    = PREC_DIR_OR_LUMPED(K,B1)
     
    lam, alpha, numbOfIter = pcgp(F,d, CP.G, e, Prec,eps0,maxIt,True,False)        
#    
    uu = []
    cnt = 0
    print('size(lam):',lam.shape)
    print('size(B):',B1[0][0].shape)
    print('type(B):',type(B1[0][0]))
    delta = 0.0
    norm_f = 0.0
    for i in range(len(K)):
        uu.append([])
        for j in range(len(K[i])):
            ind = np.arange(0,R[i][j].shape[1]) + cnt
            f_BtLam_i_j = f[i][j]-sparse.csc_matrix.dot(B1[i][j].transpose(),lam)
            KplusBtLam_i_j=Kplus[i][j]*f_BtLam_i_j
            R_alpha_i_j = sparse.csc_matrix.dot(R[i][j],alpha[ind])
            uu[i].append(KplusBtLam_i_j+R_alpha_i_j)
            cnt += R[i][j].shape[1]
            delta += np.linalg.norm(sparse.csc_matrix.dot(K[i][j],uu[i][j])-f_BtLam_i_j)                
            norm_f += np.linalg.norm(f[i][j])
        
#    u = Kplus_sub*f_m_BtLam + Roperator.mult(alpha)      
#     
#    delta = np.linalg.norm(Koperator.mult(u)-f_m_BtLam)
#    normf = np.linalg.norm(f)
    print('||Ku-f+BtLam||/||f||= %3.5e'% (np.sqrt(delta)/norm_f))
    return uu,lam
############################################################################### 
     
 
###############################################################################  
def cg(A,b,x0,R,eps0,maxIt,Prec,disp=False):
     
    n = b.shape[0]     
    x     = x0

 
    g               = A.mult(x)-b     
    Mg              = Prec*g  
    gtMg0           = np.dot(g,Mg) 
    sqrt_gtMg0      = np.sqrt(gtMg0)
    gtMg            = gtMg0
    vec_normed_g    = np.zeros(n)
 
    
    w     = Mg.copy()
     
    for i in range(n):
        
        Aw          = A.mult(w)
        rho         = -np.dot(g,Mg)/np.dot(w,Aw)
        x           += w * rho                   
        
        g           += Aw * rho
        Mg          = Prec*g      
        
        gtMgprev    = gtMg
        gtMg        = np.dot(g,Mg)
        gamma       = gtMg/gtMgprev
        w           = Mg + w * gamma 
        
        sqrt_gtMg   = np.sqrt(gtMg)
        normed_gi   = sqrt_gtMg/sqrt_gtMg0   

        vec_normed_g[i]    = normed_gi
        
        #print('....',np.log10(normed_gi/vec_normed_g[:i+1].min()),end=' ')
        if np.log10(normed_gi/vec_normed_g[:i+1].min()) > 2:
            print('... stagnate',end='')
            break        
        
        
        if disp:
            print('  ___________ inner : %d ||g||: ' % (i+2,sqrt_gtMg/sqrt_gtMg0))
        if i>maxIt or normed_gi<eps0:
            break
        
    
    numb_it = i     
    return x, sqrt_gtMg/sqrt_gtMg0, numb_it
############################################################################### 