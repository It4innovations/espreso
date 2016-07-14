#config
import multiprocessing


eps_dual_feti   = 1e-4
maxIt_dual_feti = 500


iterative_Kplus                     = False     # True, False
eps_iter_Kplus                      = 1e-15


Dirichlet_from_espreso              = False
precondDualSystem                   = 'dirichlet' # 'none', 'lumped', 'dirichlet'

single_precision                    = True     # True, False
methodToImproveSolByIterMethod      = 'pcg_x'   # 'cg_x', 'pcg_x',       
precondFrom_Areg_orA                = False
precondPrimalSystem                 = 'diag'    # 'diag', 'LU_SP', 'ILU,  'none'

#       WARNING:  ILU (incomplete LU) can generate 
#       unsymmetric preconditioner matrix!!!

# if 'cg_x' or is set
mult_Areg_or_A_RRt                  = True
flag_multiprocessing                = False



#pool = multiprocessing.Pool()




###############################################################################
###############################################################################
###############################################################################
###############################################################################
if single_precision==True:
    iterative_Kplus==True    
if precondPrimalSystem ==  'LU_SP':
    single_precision = True

#if methodToImproveSolByIterMethod = 'pcg_x':
#    precondFrom_Areg_orA = False
    
    
    
# methodToImproveSolByIterMethod      = 'cg_x'
# precondFrom_A_or_Areg               = True
# mult_Areg_or_A_RRt                  = False
#
#
#
#
    
#'cg_dx''cg_dx' 
