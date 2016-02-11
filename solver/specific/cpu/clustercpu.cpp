
#include "clustercpu.h"

void ClusterCPU::Create_SC_perDomain(bool USE_FLOAT) {

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : using Pardiso SC : ";

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		if (cluster_global_index == 1) cout << "."; // << i ;

		SparseSolverCPU tmpsps;
		if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
		tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

		if (USE_FLOAT){
			domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
			domains[i].B1Kplus.USE_FLOAT = true;
		}


	}


	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1t_comp_dom.Clear();

	if (cluster_global_index == 1)
		cout << endl;

}
