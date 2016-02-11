
#include "clusteracc.h"



void ClusterBase::Create_SC_perDomain(bool USE_FLOAT) {

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : using Pardiso SC : ";

	this->NUM_MICS = 2;
	#ifdef MIC

		// compute sizes of data to be offloaded to MIC
		eslocal maxDevNumber = this->NUM_MICS;
		if (this->NUM_MICS == 0) {
			maxDevNumber = 1;
		}
		eslocal matrixPerPack = domains.size() / maxDevNumber;
		eslocal offset = 0;
		bool symmetric = true;
		this->B1KplusPacks.resize( maxDevNumber );
		eslocal * dom2dev = new eslocal[ domains.size() ];
		eslocal * offsets = new eslocal[maxDevNumber];

		for ( eslocal i = 0; i < maxDevNumber; i++ ) {
			if ( i == maxDevNumber - 1 ) {
				matrixPerPack += domains.size() % maxDevNumber;
			}

			long dataSize = 0;
			offsets[i] = offset;

			for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
				if (!symmetric) {
					dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
				} else {
					// isPacked => is symmetric
					dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
						( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
				}
				dom2dev[ j ] = i;
			}

			this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );

			for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
				this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
					domains[j].B1t_comp_dom.cols,  symmetric );
			}
			offset += matrixPerPack;
		}
	//	tbb::mutex m;
	#endif




	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		if (cluster_global_index == 1) cout << "."; // << i ;

		SparseSolverCPU tmpsps;
		if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
		tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

		if (USE_FLOAT){
			domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
			domains[i].B1Kplus.USE_FLOAT = true;
		}

//		SparseSolverCPU tmpsps2;
//		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
//		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );

#ifdef CUDA
		//domains[i].B1Kplus.RemoveLowerDense();
		eslocal status;
		status = domains[i].B1Kplus.CopyToCUDA_Dev();
		//domains[i].B1Kplus.CopyToCUDA_Dev_fl();
		if (status == 0)
			domains[i].isOnACC = 1;
		else
			domains[i].isOnACC = 0;
#endif

#ifdef MIC
	this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
	domains[i].B1Kplus.Clear();
	//domains[i].B1t_comp_dom.Clear();
	//if (numDevices > 0) {
	//	domains[i].B1Kplus.CopyToMIC_Dev();
	//}
#endif

#ifdef CUDA
		if ( USE_KINV == 1 ) {
			cilk_for (eslocal d = 0; d < domains.size(); d++) {
				cudaError_t status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff, domains[d].B1_comp_dom.rows * sizeof(double));
				if (status != cudaSuccess)
					printf("Error allocating pinned host memory \n");

				//status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff_fl, domains[d].B1_comp_dom.rows * sizeof(float));
				//if (status != cudaSuccess)
				//	printf("Error allocating pinned host memory \n");
			}
		}
#endif

	}


#ifdef MIC
	delete [] dom2dev;
	delete [] offsets;
	if (this->NUM_MICS == 0) {
		this->B1KplusPacks[0].AllocateVectors( );
	}
	for (eslocal i = 0; i < this->NUM_MICS ; i++ ) {
		this->B1KplusPacks[i].AllocateVectors( );
		this->B1KplusPacks[i].SetDevice( i );
		this->B1KplusPacks[i].CopyToMIC();
	}

#endif

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1t_comp_dom.Clear();

	if (cluster_global_index == 1)
		cout << endl;

}
