import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import config_espreso_python as conf
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
class KPLUS:     
    # 
    def __init__(self, A, Areg, R,diagR): 
        
        self.Areg               = Areg    
        self.A                  = A.copy()
        self.iAreg              = np.ones(self.A.shape[0])
        if conf.iterative_Kplus:      
            self.R                  = R            
            #
            if conf.precondPrimalSystem == 'diag':
                if conf.precondFrom_Areg_orA:
                    self.iAreg = self.Areg.diagonal()**-1
                else:      
                    if conf.mult_Areg_or_A_RRt:
                        self.iAreg = self.A.diagonal()**-1
                    else:
                        self.iAreg = (diagR+self.A.diagonal())**-1 
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
            self.iAreg = spla.splu(self.Areg)


    def solve(self,b): 
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
                        conf.eps_iter_Kplus, maxIt, False)
            else:
                print('chose from: cg_x, cg_dx or pcg_x')
                    
            print(numbOfIter,',',end='')   
            print('%d,' % (numbOfIter),end='')
            return x            
        else:
            return self.iAreg.solve(b)  
###############################################################################       
class KPLUS_HFETI:
    def __init__(self, B0,G0,Kplus_sub,R,iF0,iS0):
        self.B0 = B0
        self.G0 = G0
        self.nB0 = B0.shape[0]
        self.Kplus_sub = Kplus_sub
        self.R = R
        self.iF0 = iF0
        self.iS0 = iS0
        
        self.B_Kplus_sub = np.zeros(B0.shape)
        B0_array = B0.toarray()
        for i in range(B0.shape[0]):
            self.B_Kplus_sub[i,:] = self.Kplus_sub.solve(B0_array[i,:])
        
    def solve(self,bc):
        b = bc[:-self.nB0:]
        c = bc[-self.nB0::]      
        Kplus_b = self.Kplus_sub.solve(b)
        d0 = sparse.csc_matrix.dot(self.B0,Kplus_b)-c 
        e0 = -sparse.csc_matrix.dot(self.R.transpose(),b)   
        G0tiF0d0 = sparse.csc_matrix.dot(self.G0.transpose(),self.iF0.solve(d0))
        beta = self.iS0.solve(G0tiF0d0-e0) 
        mu = self.iF0.solve(d0 - sparse.csc_matrix.dot(self.G0,beta))    
        x = Kplus_b-np.dot(self.B_Kplus_sub.transpose(),mu) + \
           sparse.csc_matrix.dot(self.R,beta)  
#       x = Kplus_b-self.Kplus_sub.solve(sparse.csc_matrix.dot(self.B0.transpose(),mu)) + \
#            sparse.csc_matrix.dot(self.R,beta)  
     
        return np.concatenate((x,mu))        
###############################################################################      
class FETIOPERATOR: 
    def __init__(self,Kplus_sub,B):
        self.Kplus_sub = Kplus_sub
        self.B    = B        
    def mult(self,x):
        y = sparse.csc_matrix.dot(self.B.transpose(),x)
        y = self.Kplus_sub.solve(y) 
        return sparse.csc_matrix.dot(self.B,y)
###############################################################################
class PROJ:
    def __init__(self,G,iGtG):
        self.G = G
        self.iGtG = iGtG
    def mult(self,x):
        y = x - sparse.csc_matrix.dot(self.G,\
                (self.iGtG.solve(sparse.csc_matrix.dot(self.G.transpose(),x))))
        return y        
###############################################################################
class PREC_DIR_OR_LUMPED:
    def __init__(self,K,B):
        self.K = K
        self.B = B
    def mult(self,x):
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
    tmp = np.loadtxt(pathToFile)
    I = tmp[1::,0]-1;    J = tmp[1::,1]-1;    V = tmp[1::,2]
    n = np.int32(tmp[0,0])   
    m = np.int32(tmp[0,1])
    print(str0,i,j)
    if (makeSymmetric):
        logInd = J>I; 
        I = np.concatenate((I,J[logInd]))
        J = np.concatenate((J,I[logInd]))
        V = np.concatenate((V,V[logInd]))    
    
    if (makeSparse):
        tmp = sparse.coo_matrix((V,(I,J)),shape=(n,m)).tocsc()
    else:
        if (m==1):
            tmp = V
        else:                
            tmp = sparse.coo_matrix((V,(I,J)),shape=(n,m)).toarray()
        
    return tmp
###############################################################################
def load_vector(path,str0,i,j):
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile)
    return tmp
###############################################################################  
def pcgp(F, d, G, e, Prec, eps0, maxIt,disp):
  
    
    GtG     = sparse.csc_matrix.dot(G.transpose(),G).tocsc()
    iGtG    = spla.splu(GtG)
    Proj    = PROJ(G,iGtG)
    lamIm   = sparse.csc_matrix.dot(G,iGtG.solve(e))
    nDual   = lamIm.shape[0]
    lam     = lamIm.copy()
    g       = F.mult(lam)-d     
    Pg      = Proj.mult(g) 
    MPg     = Prec.mult(Pg)
    PMPg    = Proj.mult(MPg)
    
    sqrt_gtPMPg0 = np.sqrt(np.dot(g,PMPg))
    
    if (np.dot(g,PMPg)<0): 
        raise SystemExit("Problem, precond. M is unsymmetric. Change it.") 
    sqrt_gtPg0 = np.sqrt(np.dot(g,Pg))
    
    vec_normed_g   = np.zeros(nDual)
    
#################################################    
#    PMP = np.zeros((nDual,nDual))
#    M = np.zeros((nDual,nDual))
#    e1  = np.zeros(nDual)
#    for j in range(nDual):
#        M[:,j] = Prec.mult(e1)
#        PMP[:,j]  =  Proj.mult(Prec.mult(Proj.mult(e1)))
#        if j<nDual-1:
#            e1[j]=0
#            e1[j+1]=1   
#    np.savetxt('M',M)
#    np.savetxt('PMP',PMP)
#####################################################
    
    
    w     = PMPg.copy()
    strFormat= '%3.5f'
    
    if disp:
#        print('sqrt_gtPg0: ',   sqrt_gtPg0)
#        print('sqrt_gtPMPg0: ', sqrt_gtPMPg0)
#        print('  i: ',1, '||g||: ',1)
        print('sqrt_gtPg0: %3.5e' %   (sqrt_gtPg0))
        print('sqrt_gtPMPg0:  %3.5e' % sqrt_gtPMPg0)
#        print('  i:  %d, ||g||: ',1)
#       kkk = 'i: %d, |g|: '+strFormat
        print('i: %d, |g|: %3.5f' %  (0,1))        
        
    for i in range(nDual):
        
        Fw          = F.mult(w)
        rho         = -np.dot(g,PMPg)/np.dot(w,Fw)
        lam         += w * rho        
        gprev       = g.copy()                  
        
        g           += Fw * rho
        Pg          = Proj.mult(g)
        MPg         = Prec.mult(Pg)           
        PMPgprev    = PMPg.copy()      
    

        PMPg        = Proj.mult(MPg)
        gtPMPg      = np.dot(g,PMPg)
        gamma       = gtPMPg/np.dot(gprev,PMPgprev)
        w           = PMPg + w * gamma 
        
#am_dbg        print('gtPMPg', gtPMPg)
                
        if (np.dot(g,PMPg)<0): 
            raise SystemExit("Problem, precond. M is unsymmetric. Change it.") 
          
        sqrt_gtPMPg = np.sqrt(gtPMPg)
        normed_gi   = sqrt_gtPMPg/sqrt_gtPg0
        vec_normed_g[i]    = normed_gi
        
        #print('....',np.log10(normed_gi/vec_normed_g[:i+1].min()),end=' ')
        if np.log10(normed_gi/vec_normed_g[:i+1].min()) > 2:
            print('... stagnate',end='')
            break
        
        if disp:
#            print('  i: ',i+2, '||g||: ',normed_gi)#,\
#                    #'log(||g||)/log(it): ', np.log2(normed_gi)/(i+2))
            print('i: %d, |g|: %3.5f, log(|g|)/log(it): %3.5f' % \
                                (i+2,normed_gi,np.log2(normed_gi)/(i+2)))
        if normed_gi<eps0:
            break
        if i==maxIt:
            print('PCPG does not converge within maxNumbIter (',
                       conf.maxIt_dual_feti,').')
            break
    alpha = iGtG.solve(sparse.csc_matrix.dot(G.transpose(),d-F.mult(lam)))
    numbOfIter = i    
    return lam, alpha, numbOfIter
###############################################################################  
def cg(A,b,x0,R,eps0,maxIt,Prec,disp=False):
     
    n = b.shape[0]     
    x     = x0

 
    g               = A.mult(x)-b     
    Mg              = Prec.mult(g)  
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
        Mg          = Prec.mult(g)      
        
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
def feti(K,Kreg,f,B,R,diagR,weight):
        
    weight = 1
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
              
     
    Kplus_sub=  KPLUS(K,Kreg,R,diagR)
    
      
    F       = FETIOPERATOR(Kplus_sub,B)
    G       = -sparse.csc_matrix.dot(B,R)
    e       = -sparse.csc_matrix.dot(R.transpose(),f) 
    d       = sparse.csc_matrix.dot(B,Kplus_sub.solve(f))    
    Prec    = PREC_DIR_OR_LUMPED(K,B)
    
    lam, alpha, numbOfIter = pcgp(F,d, G, e, Prec,eps0,maxIt,True)        
    
    f_m_BtLam = f - sparse.csc_matrix.dot(B.transpose(),lam)    
    u = Kplus_sub.solve(f_m_BtLam) + \
                    sparse.csc_matrix.dot(R,alpha)      
     
    delta = np.linalg.norm(sparse.csc_matrix.dot(K,u)-f_m_BtLam)
    normf = np.linalg.norm(f)
    print('||Ku-f+BtLam||/||f||= %3.5e'% (delta/normf))
    return u,lam
###############################################################################    
def hfeti(K,Kreg,f,B0,B1,R,diagR,weight):
    
    weight = 1
    maxIt = conf.maxIt_dual_feti
    eps0  = conf.eps_dual_feti   
     
    nB0         = B0.shape[0] 
    nB1         = B1.shape[0] 
    
    Oblock   = sparse.spdiags(np.zeros(nB0),0,nB0,nB0) 
    KB0t     = sparse.hstack((K,B0.transpose()))  
    B0O      = sparse.hstack((B0,Oblock)) 
    K_hfeti  = sparse.vstack((KB0t,B0O)) 
    
    Oblock2  = sparse.spdiags(np.zeros(nB1),0,nB1,nB0)
    B_hfeti  = sparse.hstack((B1,Oblock2))
    
    G0  = -sparse.csc_matrix.dot(B0,R)  
    G0tG0 = sparse.csc_matrix.dot(G0.transpose(),G0)    
    G0tG0dense = G0tG0.todense()
    U,s,V = np.linalg.svd(G0tG0dense) 
      
    
    for i in range(1,s.shape[0]):
        rat = s[i]/s[i-1]
        if np.abs(rat)<1e-6:
            H = U[:,i::]
            break
    
# -----------------------------------------------------------------------------
    # F0 create

    Kplus_sub=  KPLUS(K,Kreg,R,diagR) # Kplus_sub =  spla.splu(Kreg) 
    B0tdenseKplus_sub = B0.copy().toarray()
    
    for i in range(B0tdenseKplus_sub.shape[0]):
        B0tdenseKplus_sub[i,:] = Kplus_sub.solve(B0tdenseKplus_sub[i,:]) 
# -----------------------------------------------------------------------------
    B0dense = B0.todense()  
    iF0 = DENS_SOLVE(np.dot(B0tdenseKplus_sub,B0dense.transpose()))
   
    G0dense = G0.todense()    
    S0  =  np.dot(G0dense.transpose(), iF0.solve(G0dense))  
    for i in range(1,7):
        S0[-i::,-i::]+=S0[-1,-1]; 
    iS0 = DENS_SOLVE(S0)
    
    #iS0 = DENS_SOLVE(S0 + np.dot(H,H.transpose()))
# -----------------------------------------------------------------------------
    R_hfeti0 = sparse.csc_matrix.dot(R,H) 
    mR_hfeti0=R_hfeti0.shape[1]
    Oblock3  = sparse.spdiags(np.zeros(mR_hfeti0),0,nB0,mR_hfeti0) 
    R_hfeti  = sparse.vstack((R_hfeti0 ,Oblock3))  
    f_hfeti  = np.concatenate((f,np.zeros(nB0)))  
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    G_hfeti     = -sparse.csc_matrix.dot(B_hfeti,R_hfeti) 
    e_hfeti     = -sparse.csc_matrix.dot(R_hfeti.transpose(),f_hfeti)
    
    Kplus_hfeti = KPLUS_HFETI( B0,G0,Kplus_sub,R,iF0,iS0)
        
    
                    
    d_hfeti = sparse.csc_matrix.dot(B_hfeti,Kplus_hfeti.solve(f_hfeti))    
    
    F_hfeti = FETIOPERATOR(Kplus_hfeti,B_hfeti)
    Prec_hfeti = PREC_DIR_OR_LUMPED(K_hfeti,B_hfeti)
    
    lam_hfeti, alpha_hfeti, numbOfIter = pcgp(F_hfeti,d_hfeti, \
                        G_hfeti, e_hfeti, Prec_hfeti,eps0,maxIt,True)        
    f_m_BtLam_hfeti = f_hfeti - \
            sparse.csc_matrix.dot(B_hfeti.transpose(),lam_hfeti)
    
    u_hfeti = Kplus_hfeti.solve(f_m_BtLam_hfeti) + \
                    sparse.csc_matrix.dot(R_hfeti,alpha_hfeti)                     
                    
                    
  
                    
    
    u = u_hfeti[:-nB0]
    lam = np.concatenate((u_hfeti[-nB0:],lam_hfeti))
    
    return u, lam
     
 
 