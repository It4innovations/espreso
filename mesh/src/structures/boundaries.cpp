#include "boundaries.h"

Boundaries::Boundaries(const Mesh &mesh):
	_boundaries(mesh.coordinates().size() + mesh.coordinates().getOffset()), _mesh(mesh)
{
	const std::vector<idx_t> &parts = mesh.getPartition();
	const std::vector<Element*> &elements = mesh.getElements();

	for (size_t p = 0; p + 1 < parts.size(); p++) {
		const std::vector<idx_t> &l2g = mesh.coordinates().localToGlobal(p);
		for (idx_t e = parts[p]; e < parts[p + 1]; e++) {
			for (size_t n = 0; n < elements[e]->size(); n++) {
				_boundaries[l2g[elements[e]->node(n)]].insert(p);
			}
		}
	}
}

void Boundaries::create_B1_l(	std::vector < SparseIJVMatrix >      & B1_local, 
								std::vector < SparseIJVMatrix >      & B0_local,
								std::vector < std::vector <int> >    & l2g_vec,
								std::vector < std::vector <int> >    & lambda_map_sub_clst,
								std::vector < std::vector <int> >    & lambda_map_sub_B1,
								std::vector < std::vector <int> >    & lambda_map_sub_B0, 
								std::vector < std::vector <double> > & B1_l_duplicity,
								std::map < int, double >             & dirichlet_x,
								std::map < int, double >             & dirichlet_y,
								std::map < int, double >             & dirichlet_z,
								const int domains_num)
{

	l2g_vec.resize(domains_num);
	
	std::vector < SparseDOKMatrix > B1_loc;
	std::vector < SparseDOKMatrix > B0_loc;

	lambda_map_sub_B1.resize (domains_num); 
	lambda_map_sub_B0.resize (domains_num); 
	B1_l_duplicity.resize    (domains_num);

	int local_prim_numbering_offset = 0; 
	std::vector <int> local_prim_numbering ( domains_num, local_prim_numbering_offset);

	for (int d = 0; d < domains_num; d++) {
		int dimension = _mesh.getPartNodesCount(d) * Point::size();
		
		B1_loc.push_back( SparseDOKMatrix (1, dimension) );
		B0_loc.push_back( SparseDOKMatrix (1, dimension) );
		
		B1_local.  push_back( SparseIJVMatrix (0,0) );
		B0_local.  push_back( SparseIJVMatrix (0,0) );
	}

	std::set<int>::const_iterator it;
	std::set<int>::const_iterator it1;
	std::set<int>::const_iterator it2;
	std::map<int,double>::const_iterator itm;

	int offset = _mesh.coordinates().getOffset();
	int lambda_count_B1 = 0; 
	int lambda_count_B0 = 0; 

	for (size_t i = offset; i < _boundaries.size(); i++) {
		
		//std::vector < int > tmp_v;
		//for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
		//	tmp_v.push_back(*it);
		//}

		//if ( _boundaries[i].size() == 1 && i < 37 ) {
		//	for (int d_i = 0; d_i < 3; d_i++) {
		//		B0_loc[tmp_v[0]](lambda_count_B0, local_prim_numbering[tmp_v[0]] + d_i) =  1.0; // 3*i + d_i
		//		lambda_count_B0++;
		//	}
		//}

		//if ( _boundaries[i].size() == 2 && i < 37 ) {
		//	for (int d_i = 0; d_i < 3; d_i++) {
		//		B0_loc[tmp_v[0]](lambda_count_B0, local_prim_numbering[tmp_v[0]] + d_i) =  1.0; // 3*i + d_i
		//		lambda_count_B0++;
		//		B0_loc[tmp_v[1]](lambda_count_B0, local_prim_numbering[tmp_v[1]] + d_i) =  1.0; // 3*i + d_i
		//		lambda_count_B0++;
		//	}
		//}

		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			if ( (itm = dirichlet_x.find(i)) != dirichlet_x.end() ) {
				B1_loc           [*it](lambda_count_B1, local_prim_numbering[*it] + 0) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < int > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( (itm = dirichlet_y.find(i)) != dirichlet_y.end() ) {
				B1_loc           [*it](lambda_count_B1, local_prim_numbering[*it] + 1) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < int > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( (itm = dirichlet_z.find(i)) != dirichlet_z.end() ) {
				B1_loc           [*it](lambda_count_B1, local_prim_numbering[*it] + 2) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < int > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
		}
		
		if ( _boundaries[i].size() > 1 ) {

			// with duplicity 
			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
					for (int d_i = 0; d_i < 3; d_i++) {
						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0; 
						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0; 

						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
						lambda_map_sub_clst.push_back( std::vector < int > (1,lambda_count_B1) );

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
						for (int d_i = 0; d_i < 3; d_i++) {
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
			l2g_vec[*it].push_back( i );
		}
	
	}


	for (int d = 0; d < domains_num; d++) {
		B1_loc[d].resize(lambda_count_B1, local_prim_numbering[d]);
		B1_local[d] = SparseIJVMatrix( B1_loc[d] );

		B0_loc[d].resize(lambda_count_B0, local_prim_numbering[d]);
		B0_local[d] = SparseIJVMatrix( B0_loc[d] );
	}

}

std::ostream& operator<<(std::ostream& os, const Boundaries &b)
{
	std::set<int>::const_iterator it;
	int offset = b._mesh.coordinates().getOffset();

	for (size_t i = offset; i < b._boundaries.size(); i++) {
		os << b._mesh.coordinates()[i] << ": ";
		for (it = b._boundaries[i].begin(); it != b._boundaries[i].end(); ++it) {
			os << *it << " ";
		}
		os << "\n";
	}
	return os;
}



