
#include "boundaries.h"

namespace mesh {

template<typename T>
void Boundaries::create_B1_l(	std::vector < SparseIJVMatrix >         & B1_local,
								std::vector < SparseIJVMatrix >         & B0_local,
								std::vector < std::vector <T> >    & l2g_vec,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B0,
								std::vector < std::vector <double> >    & B1_l_duplicity,
								const eslocal domains_num)
{

	const std::map<eslocal, double> &dirichlet_x = _mesh.coordinates().property(CP::DIRICHLET_X).values();
	const std::map<eslocal, double> &dirichlet_y = _mesh.coordinates().property(CP::DIRICHLET_Y).values();
	const std::map<eslocal, double> &dirichlet_z = _mesh.coordinates().property(CP::DIRICHLET_Z).values();

	l2g_vec.resize(domains_num);

	std::vector < SparseDOKMatrix > B1_loc(domains_num);
	std::vector < SparseDOKMatrix > B0_loc(domains_num);

	lambda_map_sub_B1.resize(domains_num);
	lambda_map_sub_B0.resize(domains_num);
	B1_l_duplicity.resize(domains_num);
	std::vector<eslocal> local_prim_numbering(domains_num, 0);
	B1_local.resize(domains_num);
	B0_local.resize(domains_num);

	std::set<eslocal>::const_iterator it;
	std::set<eslocal>::const_iterator it1;
	std::set<eslocal>::const_iterator it2;

	eslocal lambda_count_B1 = 0;
	eslocal lambda_count_B0 = 0;

	for (T i = 0; i < _boundaries.size(); i++) {
		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			if ( dirichlet_x.find(index(i)) != dirichlet_x.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 0) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector <eslocal> (1, lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( dirichlet_y.find(index(i)) != dirichlet_y.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 1) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( dirichlet_z.find(index(i)) != dirichlet_z.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 2) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < eslocal > (1, lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
		}

		if ( _boundaries[i].size() > 1 ) {

			// with duplicity
			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
					for (eslocal d_i = 0; d_i < 3; d_i++) {
						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0;
						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0;

						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
						lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );

						B1_l_duplicity[*it1].push_back( 1.0 / (double)_boundaries[i].size() );
						B1_l_duplicity[*it2].push_back( 1.0 / (double)_boundaries[i].size() );

						lambda_count_B1++;
					}
				}
			}


			// no duplicity
			bool is_corner = true;
			if ( is_corner ) {
				for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
					if (it1 != _boundaries[i].begin()) {
						for (eslocal d_i = 0; d_i < 3; d_i++) {
							B0_loc[*it2](lambda_count_B0, local_prim_numbering[*it2] + d_i) =  1.0;
							B0_loc[*it1](lambda_count_B0, local_prim_numbering[*it1] + d_i) = -1.0;
							lambda_map_sub_B0[*it2].push_back(lambda_count_B0);
							lambda_map_sub_B0[*it1].push_back(lambda_count_B0);
							lambda_count_B0++;
						}
					}
					it2 = it1;
				}
			}
		}
		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			local_prim_numbering[*it] += 3;
			l2g_vec[*it].push_back(index(i));
		}
	}

	for (eslocal d = 0; d < domains_num; d++) {
		B1_loc[d].resize(lambda_count_B1, local_prim_numbering[d]);
		B1_local[d] = B1_loc[d];

		B0_loc[d].resize(lambda_count_B0, local_prim_numbering[d]);
		B0_local[d] = B0_loc[d];
	}
}

}



