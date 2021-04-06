
#include "geometry.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactstore.h"
#include "math/matrix.dense.h"

#include <vector>

namespace espreso {
namespace geometry {

void computeBoundaryRegionsArea()
{
	std::vector<double> area(info::mesh->boundaryRegions.size()), garea(info::mesh->boundaryRegions.size());

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			BoundaryRegionStore *store = info::mesh->boundaryRegions[r];

			double A = 0;
			auto nodes = store->procNodes->cbegin();
			const auto &epointers = store->epointers->datatarray();
			const auto &coordinates = info::mesh->nodes->coordinates->datatarray();
			for (size_t e = 0; e < store->procNodes->structures(); ++e, ++nodes) {

				MatrixDense coords(nodes->size(), 3), dND(1, 3);

				const std::vector<MatrixDense> &dN = *epointers[e]->dN;
				const std::vector<double> &weighFactor = *epointers[e]->weighFactor;

				for (size_t n = 0; n < nodes->size(); ++n) {
					coords(n, 0) = coordinates[nodes->at(n)].x;
					coords(n, 1) = coordinates[nodes->at(n)].y;
					coords(n, 2) = coordinates[nodes->at(n)].z;
				}

				if (store->dimension == 1) {
					for (size_t gp = 0; gp < dN.size(); gp++) {
						dND.multiply(dN[gp], coords);
						A += dND.norm() * weighFactor[gp];
					}
				}
				if (store->dimension == 2) {
					for (size_t gp = 0; gp < dN.size(); gp++) {
						dND.multiply(dN[gp], coords);
						Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
						Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
						Point va = Point::cross(v1, v2);
						A += va.norm() * weighFactor[gp];
					}
				}
			}
			area[r] = A;
		}
	}

	Communication::allReduce(area.data(), garea.data(), area.size(), MPI_DOUBLE, MPI_SUM);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			info::mesh->boundaryRegions[r]->area = garea[r];
		}
	}
}

void assembleMortarInterface(std::vector<ijv> &B)
{
	// Popp: Mortar Methods for Computational Contact Mechanics and General Interface Problems
	// dissertation at Technische universitat Munchen
	// https://mediatum.ub.tum.de/doc/1109994/1109994.pdf

//	int printrank = -1;
//	std::vector<int>    slJs;   slJs.reserve( 4);
//	int maJoffset = 0;
//	std::vector<int>    maJs;   maJs.reserve(40);
//	// return values
//	std::vector<ijv> Be;
//	//                segment = [ surface element offset; number of slaves; neighbor offset; slave offset ]
//	// segment = [ surface element offset; number of master elements; neighbor offset; master offset ]
//	auto triangles = _mesh->contacts->intersections->begin();
//	std::vector<double> slCover(_mesh->contacts->contactPairs->structures(), 0.0);
//	esint segmentIndex = 0;
//	std::vector<double> me;   me.reserve(9*9);
//	std::vector<double> de;   de.reserve(9*9);
//	size_t slNodesSize = 0, maNodesSize = 0, tmpCnt;
//
//	struct dmAssembleDataForGPItem {
//		double wJtJs;
//		int maJoffset;
//		std::vector<double> slN;
//		std::vector<double> maN;
//	};
//	std::vector<dmAssembleDataForGPItem> dmAssembleDataForGP; dmAssembleDataForGP.reserve(40*1);  // CORRECT WITH ONDRA
//	if (info::mpi::rank == printrank) std::cout << "intersections: " << _mesh->contacts->intersections->structures() << std::endl;
//	if (info::mpi::rank == printrank) std::cout << "triangles    : " << _mesh->contacts->intersections->datatarray().size() << std::endl;
//
//	Element *trElement = _mesh->edata[(int)Element::CODE::TRIANGLE3];
//	size_t trQps = trElement->weighFactor->size();
//	MatrixDense trCoordsMat(3,2), tmpMatrixDense;
//	MatrixDense trDetJ(trQps,1), slDetJ(trQps,1), maDetJ(trQps,1), slN, maN;
//	if (info::mpi::rank == printrank) std::cout << "trQps: " << trQps << std::endl;
//	MatrixDense slTrNodRefPoints(3,2), slTrQpRefPoints(trQps,2), maTrNodRefPoints(3,2), maTrQpRefPoints(trQps,2);
//	for (auto segment= _mesh->contacts->contactPairs->begin(); segment != _mesh->contacts->contactPairs->end(); ++segment, ++segmentIndex) {
//		slJs.resize(0);
//		maJs.resize(0);
//		if (info::mpi::rank == printrank) std::cout << "slave [";
//		std::vector<Point> slCoords;
//		auto slNodesIt = _mesh->surface->enodes->begin() + segment->at(0);
//		Element* slElement = _mesh->surface->epointers->datatarray()[segment->at(0)];
//
//		MatrixDense slCoordsMat(slNodesIt->size(),2);
//		dmAssembleDataForGP.reserve(5*segment->size());                               // CORRECT WITH ONDRA
//		dmAssembleDataForGP.resize(0);
//
//		slNodesSize = slElement->nodes;
//		for (auto n = slNodesIt->begin(); n != slNodesIt->end(); ++n) {
//			if (info::mpi::rank == printrank) std::cout << " " << *n;
//			slJs.push_back(*n);
//		}
//		de.resize(slNodesSize,slNodesSize);
//		me.resize(slNodesSize,slNodesSize);
//		if (info::mpi::rank == printrank) std::cout << " ][";
//		unsigned int tmpcount = 0;
//		for (auto n = slNodesIt->begin(); n != slNodesIt->end(); ++n, ++tmpcount) {
//			if (info::mpi::rank == printrank) std::cout << " " << _mesh->nodes->coordinates->datatarray()[*n];
//			slCoords.push_back(_mesh->nodes->coordinates->datatarray()[*n]);
//			slCoordsMat[tmpcount][0] = _mesh->nodes->coordinates->datatarray()[*n].x;
//			slCoordsMat[tmpcount][1] = _mesh->nodes->coordinates->datatarray()[*n].y;
//		}
//		if (info::mpi::rank == printrank) std::cout << " ]\n";
//
//		maJoffset = 0;
//		for (esint m = 0; m < segment->at(1); ++m, ++triangles) {
//			if (_mesh->contacts->gneighbors[segment->at(2 + 2 * m)] < info::mpi::rank) {
//				continue;
//			}
//			if (info::mpi::rank == printrank) std::cout << "  master [";
//			auto maNodesIt = _mesh->contacts->gnenodes[segment->at(2 + 2 * m)]->begin() + segment->at(2 + 2 * m + 1);
//			Element* maElement = _mesh->contacts->gnepointers[segment->at(2 + 2 * m)]->datatarray()[segment->at(2 + 2 * m + 1)];
//			MatrixDense maCoordsMat(maNodesIt->size(),2);
//			maNodesSize = maElement->nodes;
//			auto maCoords = _mesh->contacts->gnecoords[segment->at(2 + 2 * m)]->begin() + segment->at(2 + 2 * m + 1);
//			for (auto n = maNodesIt->begin(); n != maNodesIt->end(); ++n) {
//				if (info::mpi::rank == printrank) std::cout << " " << *n;
//				maJs.push_back(*n);
//			}
//			if (info::mpi::rank == printrank) std::cout << " ][";
//			tmpcount = 0;
//			for (auto n = maCoords->begin(); n != maCoords->end(); ++n, ++tmpcount) {
//				if (info::mpi::rank == printrank) std::cout << " " << *n;
//				maCoordsMat[tmpcount][0] = (*n).x;
//				maCoordsMat[tmpcount][1] = (*n).y;
//			}
//			if (info::mpi::rank == printrank) std::cout << " ]\n";
//
//			if (info::mpi::rank == printrank) std::cout << "    intersetions: " << std::endl;
//			for (size_t t = 0; t < triangles->size(); ++t) {
//				// compute referenceCoordinatesMS and masterCover
//				trCoordsMat[0][0] = triangles->at(t).p[0].x;  trCoordsMat[0][1] = triangles->at(t).p[0].y;
//				trCoordsMat[1][0] = triangles->at(t).p[1].x;  trCoordsMat[1][1] = triangles->at(t).p[1].y;
//				trCoordsMat[2][0] = triangles->at(t).p[2].x;  trCoordsMat[2][1] = triangles->at(t).p[2].y;
//
//				slElement->computeReferenceCoords(slCoordsMat, trCoordsMat, slTrNodRefPoints);
//				maElement->computeReferenceCoords(maCoordsMat, trCoordsMat, maTrNodRefPoints);
//				for (size_t qp = 0; qp < trQps; ++qp) {
//					tmpMatrixDense.multiply(trElement->N->at(qp), slTrNodRefPoints);
//					slTrQpRefPoints[qp][0] = tmpMatrixDense[0][0];   slTrQpRefPoints[qp][1] = tmpMatrixDense[0][1];
//					tmpMatrixDense.multiply(trElement->N->at(qp), maTrNodRefPoints);
//					maTrQpRefPoints[qp][0] = tmpMatrixDense[0][0];   maTrQpRefPoints[qp][1] = tmpMatrixDense[0][1];
//				}
//				slCover[segmentIndex] += Point::cross2d(Point(slTrNodRefPoints[1][0]-slTrNodRefPoints[0][0],slTrNodRefPoints[1][1]-slTrNodRefPoints[0][1],0), Point(slTrNodRefPoints[2][0]-slTrNodRefPoints[0][0],slTrNodRefPoints[2][1]-slTrNodRefPoints[0][1],0));
//				trElement->recomputeDetJ( slTrNodRefPoints, trDetJ);
//				slElement->recomputeDetJN(slCoordsMat, slDetJ, slN, slTrQpRefPoints);
//				maElement->recomputeDetJN(maCoordsMat, maDetJ, maN, maTrQpRefPoints);
//				if (info::mpi::rank == printrank) std::cout << "       " << triangles->at(t);
//				if (info::mpi::rank == printrank) std::cout << " sl[<" << slTrNodRefPoints[0][0] << "," << slTrNodRefPoints[0][1] << "><" << slTrNodRefPoints[1][0] << "," << slTrNodRefPoints[1][1] << "><" << slTrNodRefPoints[2][0] << "," << slTrNodRefPoints[2][1] << ">] ma["   << maTrNodRefPoints[0][0] << "," << maTrNodRefPoints[0][1] << "><" << maTrNodRefPoints[1][0] << "," << maTrNodRefPoints[1][1] << "><" << maTrNodRefPoints[2][0] << "," << maTrNodRefPoints[2][1] << ">]" << std::endl;
//				if (info::mpi::rank == printrank) std::cout << "            slQp[<" << slTrQpRefPoints[0][0] << "," << slTrQpRefPoints[0][1]<< ">]  maQp[<" << maTrQpRefPoints[0][0] << "," << maTrQpRefPoints[0][1]<< ">]  " << std::endl;
//				for (size_t jj = 0; jj < trElement->weighFactor->size(); ++jj) {
//					dmAssembleDataForGP.push_back(dmAssembleDataForGPItem());
//					dmAssembleDataForGP.back().wJtJs = trElement->weighFactor->at(jj) * trDetJ[0][jj] * slDetJ[0][jj];
//					dmAssembleDataForGP.back().slN.resize(slElement->nodes);  for (auto ii = 0; ii < slElement->nodes; ii++) { dmAssembleDataForGP.back().slN[ii] =  slN[ii][jj]; }
//					dmAssembleDataForGP.back().maN.resize(maElement->nodes);  for (auto ii = 0; ii < maElement->nodes; ii++) { dmAssembleDataForGP.back().maN[ii] =  maN[ii][jj]; }
//					dmAssembleDataForGP.back().maJoffset = maJoffset;
//				}
//			}
//			maJoffset += maElement->nodes;
//			if (info::mpi::rank == printrank) std::cout << "\n";
//		}
//		switch (slElement->code) {
//		case Element::CODE::TRIANGLE3: slCover[segmentIndex]/=0.5; break;
//		case Element::CODE::SQUARE4:   slCover[segmentIndex]/=4.0; break;
//		default: break;
//		}
//		if (slCover[segmentIndex] >= MIN_SLAVE_COVER_RATIO) {
//			for (size_t ij = 0; ij < slNodesSize*slNodesSize; ij++) { de[ij] = 0.0; }  // clear De
//			for (size_t ij = 0; ij < slNodesSize*maNodesSize; ij++) { me[ij] = 0.0; }  // clear Me
//			if (slCover[segmentIndex] >= TOLERATED_SLAVE_COVER_RATIO_FOR_DUAL_SHAPE_COEFFICIENTS_ON_WHOLE_ELEMENT) {
//				// compute De, Me for coefficients Ae of dual basis functions from formula (4.57)
//				MatrixDense detJ;
//				slElement->recomputeDetJ(slCoordsMat, detJ);
//				for (size_t i = 0; i < slNodesSize; i++) {          // LOCAL SHAPE FUNCTIONS - i index
//					for (size_t j = i; j < slNodesSize; j++) {      // LOCAL SHAPE FUNCTIONS - j index
//						for (size_t gp = 0; gp < slElement->N->size(); gp++) {
//							me[i*slNodesSize+j]     += slElement->N->at(i)[0][gp] * slElement->N->at(j)[0][gp] * detJ[0][gp] * slElement->weighFactor->at(gp); // me = M_{e}^T == M_{e}
//							if (i == j) {
//								de[i*slNodesSize+j] += slElement->N->at(i)[0][gp] * detJ[0][gp] * slElement->weighFactor->at(gp); // me = M_{e}^T == M_{e}
//							}
//						}
//					}
//					for (size_t j = 0; j < i; j++) {
//						me[i*slNodesSize+j] = me[j*slNodesSize+i];
//					}
//				}
//			} else {
//				// compute De, Me for coefficients Ae of dual basis functions from formula (4.60)
//				for (auto it = dmAssembleDataForGP.begin(); it != dmAssembleDataForGP.end(); it++) {
//					for (size_t i = 0; i < slNodesSize; i++) {          // LOCAL SHAPE FUNCTIONS - i index
//						for (size_t j = i; j < slNodesSize; j++) {      // LOCAL SHAPE FUNCTIONS - j index
//							me[i*slNodesSize+j] += it->wJtJs * it->slN[i] * it->slN[j];
//						}
//						de[i*slNodesSize+i] += it->wJtJs * it->slN[i];
//					}
//				}
//				for (size_t i = 0; i < slNodesSize; i++) {
//					for (size_t j = 0; j < i; j++) {
//						me[i*slNodesSize+j] = me[j*slNodesSize+i];
//					}
//				}
//			}
//			// de = me^{-1}*de == M_{e}^{-T}*D_{e}^{T} == D_{e}*M_{e}^{-1} Popp (4.57)
//			MATH::DenseMatDenseMatRowMajorSystemSolve(slNodesSize, slNodesSize, me.data(), de.data());
//			// remark: slJs.size() == slNodesSize
//			Be.reserve(slNodesSize*(slNodesSize+maJs.size()));
//			tmpCnt = 0;
//			for (size_t i = 0; i < slNodesSize; i++) {
//				for (size_t j = 0; j < slNodesSize; ++j, ++tmpCnt) {
//					Be.push_back(ijv(slJs[i], slJs[j], 0.0));
//				}
//				for (size_t j = 0; j < maJs.size(); ++j, ++tmpCnt) {
//					Be.push_back(ijv(slJs[i], -1, 0.0));
//				}
//			}
////			std::cout << "dmAssembleDataForGP" << std::endl;
//			for (auto it = dmAssembleDataForGP.begin(); it != dmAssembleDataForGP.end(); it++) {
//				MatrixDense psi;
//				psi.multiply(de.data(), slNodesSize, slNodesSize, it->slN.size(), 1, it->slN.data());
//				tmpCnt = 0;
//				for (size_t i = 0; i < slNodesSize; i++) {
//					for (size_t j = 0; j < slNodesSize; ++j, ++tmpCnt) {
//						Be[tmpCnt].v                 += it->wJtJs * it->slN[j] * psi[0][i];
//					}
//					for (size_t j = 0; j < it->maN.size(); ++j) {
//						Be[tmpCnt+j+it->maJoffset].j  = maJs[j+it->maJoffset];
//						Be[tmpCnt+j+it->maJoffset].v -= it->wJtJs * it->maN[j] * psi[0][i];
//					}
//					tmpCnt += maJs.size();
//				}
//			} // end dmAssembleDataForGP
//		} // else (slCover[segmentIndex] < MIN_SLAVE_COVER_RATIO) ---> do nothing
//		for (auto itBe = Be.begin(); itBe < Be.end(); ++itBe) {
//			if ( std::abs(itBe->v) > BE_VALUE_TRESHOLD) {
//				if (info::mpi::rank == printrank) std::cout << itBe->i << "," << itBe->j << "="  << itBe->v << "\n";
//				B.push_back(*itBe);
//			}
//		}
//		Be.clear();
//	}
}

void computeMortars()
{
//	std::vector<ijv> B;
//	if (_mesh->contacts->contactPairs != NULL) {
//		assembleMortarInterface(B);
//	}
//
//	std::sort(B.begin(), B.end());
//
//	// exchange neighbors
//	std::vector<int> neighs = _mesh->contacts->gneighbors;
//	std::vector<std::vector<int> > rNeighs(_mesh->neighbors.size());
//	if (!Communication::exchangeUnknownSize(_mesh->contacts->gneighbors, rNeighs, _mesh->neighbors)) {
//		eslog::error("ESPRESO internal error: exchange mortar neighbors.\n");
//	}
//	for (size_t r = 0; r < rNeighs.size(); ++r) {
//		neighs.insert(neighs.end(), rNeighs[r].begin(), rNeighs[r].end());
//	}
//	rNeighs.resize(_mesh->contacts->gneighbors.size());
//	if (!Communication::exchangeUnknownSize(_mesh->neighbors, rNeighs, _mesh->contacts->gneighbors)) {
//		eslog::error("ESPRESO internal error: exchange mortar geometric neighbors.\n");
//	}
//	for (size_t r = 0; r < rNeighs.size(); ++r) {
//		neighs.insert(neighs.end(), rNeighs[r].begin(), rNeighs[r].end());
//	}
//	utils::sortAndRemoveDuplicates(neighs);
//	auto myrank = std::lower_bound(neighs.begin(), neighs.end(), info::mpi::rank);
//	if (myrank != neighs.end() && *myrank == info::mpi::rank) {
//		neighs.erase(myrank);
//	}
//
//	// synchronize across neighbors
//	std::vector<std::vector<ijv> > Bs(neighs.size(), B), Br(neighs.size());
//
//	if (!Communication::exchangeUnknownSize(Bs, Br, neighs)) {
//		eslog::error("ESPRESO internal error: cannot synchronize mortars.\n");
//	}
//
//	for (size_t n = 0; n < Br.size(); ++n) {
//		B.insert(B.end(), Br[n].begin(), Br[n].end());
//	}
//
//	std::sort(B.begin(), B.end());
//	if (B.size()) {
//		size_t unique = 0;
//		for (size_t i = 1; i < B.size(); i++) {
//			if (B[unique] != B[i]) {
//				B[++unique] = B[i];
//			} else {
//				B[unique].v += B[i].v;
//			}
//		}
//		B.resize(unique + 1);
//	}
//
//	neighs.push_back(info::mpi::rank);
//	std::sort(neighs.begin(), neighs.end());
//	std::vector<esint> cIDs, crankdist;
//	std::vector<int> crankdata;
//	std::vector<std::vector<esint> > sBuffer(neighs.size()), rBuffer(neighs.size());
//	for (auto begin = B.begin(), end = begin; begin != B.end(); begin = end) {
//		while (end != B.end() && begin->i == end->i) {
//			++end;
//		}
//		auto &ids = _mesh->nodes->IDs->datatarray();
//		size_t bsize = sBuffer[0].size();
//		for (auto it = begin; it != end; ++it) {
//			auto nit = std::find(ids.begin(), ids.begin() + _mesh->nodes->uniqInfo.nhalo, it->j);
//			if (nit == ids.begin() + _mesh->nodes->uniqInfo.nhalo || *nit != it->j) {
//				nit = std::lower_bound(ids.begin() + _mesh->nodes->uniqInfo.nhalo, ids.end(), it->j);
//			}
//			if (nit != ids.end() && *nit == it->j) {
//				sBuffer[0].push_back(it->j);
//			}
//		}
//		if (bsize != sBuffer[0].size()) {
//			_mesh->contacts->B.insert(_mesh->contacts->B.end(), begin, end);
//			for (auto it = begin; it != end; ++it) {
//				cIDs.push_back(it->j);
//			}
//		}
//	}
//
//	utils::sortAndRemoveDuplicates(sBuffer[0]);
//	for (size_t n = 1; n < neighs.size(); ++n) {
//		sBuffer[n] = sBuffer[0];
//	}
//
//	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighs)) {
//		eslog::error("ESPRESO internal error: cannot exchange requested mortars IDs.\n");
//	}
//
//	utils::sortAndRemoveDuplicates(cIDs);
//
//	std::vector<std::vector<esint>::const_iterator> rPointer(rBuffer.size());
//	for (size_t r = 0; r < rBuffer.size(); r++) {
//		rPointer[r] = rBuffer[r].begin();
//	}
//	crankdist.push_back(0);
//	for (size_t n = 0; n < cIDs.size(); ++n) {
//		for (size_t r = 0; r < rBuffer.size(); r++) {
//			while (rPointer[r] != rBuffer[r].end() && *rPointer[r] < cIDs[n]) {
//				++rPointer[r];
//			}
//			if (rPointer[r] != rBuffer[r].end() && *rPointer[r] == cIDs[n]) {
//				crankdata.push_back(neighs[r]);
//				++rPointer[r];
//			}
//		}
//		crankdist.push_back(crankdata.size());
//	}
//
//	_mesh->contacts->nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, cIDs));
//	_mesh->contacts->nodes->distribution = _mesh->contacts->nodes->IDs->datatarray().distribution();
//	_mesh->contacts->nodes->ranks = new serializededata<esint, int>(crankdist, crankdata);
//
//	utils::sortAndRemoveDuplicates(crankdata);
//	_mesh->neighborsWithMe.insert(_mesh->neighborsWithMe.end(), crankdata.begin(), crankdata.end());
//	utils::sortAndRemoveDuplicates(_mesh->neighborsWithMe);
//	_mesh->neighbors = _mesh->neighborsWithMe;
//	_mesh->neighbors.erase(std::lower_bound(_mesh->neighbors.begin(), _mesh->neighbors.end(), info::mpi::rank));
//
//	Communication::serialize([&] () {
////		printf(" == %d == (%ld) \n", info::mpi::rank, B.size());
////		std::cout << _mesh->neighbors;
////		std::cout << _mesh->neighborsWithMe;
////		std::cout << *_mesh->contacts->nodes->IDs << "\n";
////		std::cout << *_mesh->contacts->nodes->ranks << "\n";
////		for (auto it = _mesh->contacts->B.begin(); it != _mesh->contacts->B.end(); ++it) {
////			printf("<%2d,%2d> = %+8.6f\n", it->i, it->j, it->v);
////		}
//	});
}

}
}
