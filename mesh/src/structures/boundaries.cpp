#include "boundaries.h"

using namespace mesh;

Boundaries::Boundaries(const Mesh &m): _boundaries(m.coordinates().size()), _mesh(m)
{
	const std::vector<esint> &parts = m.getPartition();
	const std::vector<Element*> &elements = m.getElements();
	const Coordinates &c = m.coordinates();

	for (size_t p = 0; p + 1 < parts.size(); p++) {
		for (esint e = parts[p]; e < parts[p + 1]; e++) {
			for (size_t n = 0; n < elements[e]->size(); n++) {
				_boundaries[c.clusterIndex(elements[e]->node(n), p)].insert(p);
			}
		}
	}
}

void Boundaries::create_B1_l(	std::vector < SparseIJVMatrix >      & B1_local, 
								std::vector < SparseIJVMatrix >      & B0_local,
								std::vector < std::vector <esint> >    & l2g_vec,
								std::vector < std::vector <esint> >    & lambda_map_sub_clst,
								std::vector < std::vector <esint> >    & lambda_map_sub_B1,
								std::vector < std::vector <esint> >    & lambda_map_sub_B0,
								std::vector < std::vector <double> > & B1_l_duplicity,
								std::map < esint, double >             & dirichlet_x,
								std::map < esint, double >             & dirichlet_y,
								std::map < esint, double >             & dirichlet_z,
								const esint domains_num)
{

	l2g_vec.resize(domains_num);
	
	std::vector < SparseDOKMatrix > B1_loc;
	std::vector < SparseDOKMatrix > B0_loc;

	lambda_map_sub_B1.resize (domains_num); 
	lambda_map_sub_B0.resize (domains_num); 
	B1_l_duplicity.resize    (domains_num);

	esint local_prim_numbering_offset = 0;
	std::vector <esint> local_prim_numbering ( domains_num, local_prim_numbering_offset);

	for (esint d = 0; d < domains_num; d++) {
		esint dimension = _mesh.getPartNodesCount(d) * Point::size();
		
		B1_loc.push_back( SparseDOKMatrix (1, dimension) );
		B0_loc.push_back( SparseDOKMatrix (1, dimension) );
		
		B1_local.  push_back( SparseIJVMatrix (0,0) );
		B0_local.  push_back( SparseIJVMatrix (0,0) );
	}

	std::set<esint>::const_iterator it;
	std::set<esint>::const_iterator it1;
	std::set<esint>::const_iterator it2;
	std::map<esint,double>::const_iterator itm;

	esint lambda_count_B1 = 0;
	esint lambda_count_B0 = 0;

	for (size_t i = 0; i < _boundaries.size(); i++) {
		
		//std::vector < esint > tmp_v;
		//for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
		//	tmp_v.push_back(*it);
		//}

		//if ( _boundaries[i].size() == 1 && i < 37 ) {
		//	for (esint d_i = 0; d_i < 3; d_i++) {
		//		B0_loc[tmp_v[0]](lambda_count_B0, local_prim_numbering[tmp_v[0]] + d_i) =  1.0; // 3*i + d_i
		//		lambda_count_B0++;
		//	}
		//}

		//if ( _boundaries[i].size() == 2 && i < 37 ) {
		//	for (esint d_i = 0; d_i < 3; d_i++) {
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
				lambda_map_sub_clst.push_back( std::vector < esint > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( (itm = dirichlet_y.find(i)) != dirichlet_y.end() ) {
				B1_loc           [*it](lambda_count_B1, local_prim_numbering[*it] + 1) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < esint > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( (itm = dirichlet_z.find(i)) != dirichlet_z.end() ) {
				B1_loc           [*it](lambda_count_B1, local_prim_numbering[*it] + 2) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < esint > (1,lambda_count_B1) );
				B1_l_duplicity   [*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
		}
		
		if ( _boundaries[i].size() > 1 ) {

			// with duplicity 
			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
					for (esint d_i = 0; d_i < 3; d_i++) {
						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0; 
						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0; 

						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
						lambda_map_sub_clst.push_back( std::vector < esint > (1,lambda_count_B1) );

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
						for (esint d_i = 0; d_i < 3; d_i++) {
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


	for (esint d = 0; d < domains_num; d++) {
		B1_loc[d].resize(lambda_count_B1, local_prim_numbering[d]);
		B1_local[d] = SparseIJVMatrix( B1_loc[d] );

		B0_loc[d].resize(lambda_count_B0, local_prim_numbering[d]);
		B0_local[d] = SparseIJVMatrix( B0_loc[d] );
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Boundaries &b)
{
	std::set<esint>::const_iterator it;

	for (size_t i = 0; i < b._boundaries.size(); i++) {
		os << b._mesh.coordinates()[i] << ": ";
		for (it = b._boundaries[i].begin(); it != b._boundaries[i].end(); ++it) {
			os << *it << " ";
		}
		os << "\n";
	}
	return os;
}



