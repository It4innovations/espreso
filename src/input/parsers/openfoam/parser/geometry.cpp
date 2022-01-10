
#include "geometry.h"

#include "basis/containers/tarray.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

static void ofreduce(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		*(reinterpret_cast<FoamFile*>(out) + i) ^= *(reinterpret_cast<FoamFile*>(in) + i);
	}
}

void OpenFOAMGeometry::scan()
{
	points.scan(pointsFile);
	faces.scan(facesFile);
//	owner.scan(_files.owner.file);
//	neighbour.scan(_files.neighbour.file);

	MPI_Datatype oftype;
	MPI_Type_contiguous(sizeof(FoamFile), MPI_BYTE, &oftype);
	MPI_Type_commit(&oftype);

	MPI_Op ofop;
	MPI_Op_create(ofreduce, 1, &ofop);

	FoamFile files[] = { points, faces };
	Communication::allReduce(files, NULL, 2, oftype, ofop);

	MPI_Op_free(&ofop);
	MPI_Type_free(&oftype);

	static_cast<FoamFile&>(points) = files[0];
	static_cast<FoamFile&>(faces) = files[1];
}

void OpenFOAMGeometry::parse(OrderedUniqueNodes *nodes, OrderedUniqueFaces *faces)
{
	this->points.parse(nodes->coordinates);
	this->faces.parse(faces->etype, faces->enodes);
}

//static void exchangeConnectivity(size_t nfaces, ivector<esint> &owner, ivector<esint> &neighbour)
//{
////	std::vector<size_t> offset = { owner.size(), neighbour.size() };
////	std::vector<size_t> distribution = Communication::getDistribution(nfaces);
////	Communication::exscan(offset);
////
////	std::vector<esint> sBuffer, rBuffer;
//}
//
//static void exchangeFaces(ivector<esint> &fdist, ivector<esint> &fnodes, ivector<esint> &owner, ivector<esint> &neighbour)
//{
//
//}
//
//static void buildElements(ivector<esint> &fdist, ivector<esint> &fnodes, ivector<esint> &owner, ivector<esint> &neighbour, OrderedUniqueElements &elements)
//{
//	struct __element__ {
//		int triangles = 0, squares = 0, others = 0;
//
//		void add(const esint &nodes)
//		{
//			switch (nodes) {
//			case 3: ++triangles; break;
//			case 4: ++squares; break;
//			default: ++others;
//			}
//		}
//
//		Element::CODE type()
//		{
//			if (triangles == 4 && squares == 0 && others == 0) { return Element::CODE::TETRA4; }
//			if (triangles == 4 && squares == 1 && others == 0) { return Element::CODE::PYRAMID5; }
//			if (triangles == 2 && squares == 3 && others == 0) { return Element::CODE::PRISMA6; }
//			if (triangles == 0 && squares == 6 && others == 0) { return Element::CODE::HEXA8; }
//			return Element::CODE::NOT_SUPPORTED;
//		}
//	};
//
//	auto addTetra = [&] (ivector<esint>::const_iterator obegin, ivector<esint>::const_iterator oend, ivector<esint>::const_iterator nbegin, ivector<esint>::const_iterator nend) {
//
//	};
//
//	auto addPyramid = [&] (ivector<esint>::const_iterator obegin, ivector<esint>::const_iterator oend, ivector<esint>::const_iterator nbegin, ivector<esint>::const_iterator nend) {
//
//	};
//
//	auto addPrisma = [&] (ivector<esint>::const_iterator obegin, ivector<esint>::const_iterator oend, ivector<esint>::const_iterator nbegin, ivector<esint>::const_iterator nend) {
//
//	};
//
//	auto addHexa = [&] (ivector<esint>::const_iterator obegin, ivector<esint>::const_iterator oend, ivector<esint>::const_iterator nbegin, ivector<esint>::const_iterator nend) {
//		elements.etype.push_back(Element::CODE::HEXA8);
//		elements.enodes.insert(elements.enodes.end(), 8, -1); // better solution
//		auto nodes = elements.enodes.data() + elements.enodes.size() - 8;
//		if (obegin != oend) { // there is at least one owner
//			nodes[0] = fnodes[fdist[owner[*obegin]] + 0];
//			nodes[1] = fnodes[fdist[owner[*obegin]] + 1];
//			nodes[5] = fnodes[fdist[owner[*obegin]] + 2];
//			nodes[4] = fnodes[fdist[owner[*obegin]] + 3];
//		} else {
//			nodes[0] = fnodes[fdist[neighbour[*nbegin]] + 0];
//			nodes[1] = fnodes[fdist[neighbour[*nbegin]] + 1];
//			nodes[5] = fnodes[fdist[neighbour[*nbegin]] + 2];
//			nodes[4] = fnodes[fdist[neighbour[*nbegin]] + 3];
//		}
//	};
//
//	auto addOther = [&] (ivector<esint>::const_iterator obegin, ivector<esint>::const_iterator oend, ivector<esint>::const_iterator nbegin, ivector<esint>::const_iterator nend) {
//
//	};
//
//	std::vector<esint> edistribution = tarray<esint>::distribute(info::env::OMP_NUM_THREADS, fdist.size() - 1);
//
//	ivector<esint> opermutation(owner.size()), npermutation(neighbour.size());
//	std::iota(opermutation.begin(), opermutation.end(), 0);
//	std::sort(opermutation.begin(), opermutation.end(), [&owner] (const esint &i, const esint &j) { return owner[i] < owner[j]; });
//	std::iota(npermutation.begin(), npermutation.end(), 0);
//	std::sort(npermutation.begin(), npermutation.end(), [&neighbour] (const esint &i, const esint &j) { return neighbour[i] < neighbour[j]; });
//
////	#pragma omp parallel for
//	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
//		auto obegin = std::lower_bound(opermutation.cbegin(), opermutation.cbegin(), edistribution[t], [&owner] (const esint &p, const esint &e) { return owner[p] < e; });
//		auto nbegin = std::lower_bound(npermutation.cbegin(), npermutation.cbegin(), edistribution[t], [&neighbour] (const esint &p, const esint &e) { return neighbour[p] < e; });
//		auto oit = obegin, nit = nbegin;
//
//		for (esint e = edistribution[t]; oit != opermutation.end() && nit != npermutation.end(); ++e, obegin = oit, nbegin = nit) {
//			__element__ element;
//			while (oit != opermutation.end() && owner[*oit] == e) { element.add(fdist[*oit + 1] - fdist[*oit]); ++oit; }
//			while (nit != npermutation.end() && neighbour[*nit] == e) { element.add(fdist[*nit + 1] - fdist[*nit]); ++nit; }
//			switch (element.type()) {
//			case Element::CODE::TETRA4: addTetra(obegin, oit, nbegin, nit); break;
//			case Element::CODE::PYRAMID5: addPyramid(obegin, oit, nbegin, nit); break;
//			case Element::CODE::PRISMA6: addPrisma(obegin, oit, nbegin, nit); break;
//			case Element::CODE::HEXA8: addHexa(obegin, oit, nbegin, nit); break;
//			case Element::CODE::NOT_SUPPORTED: addOther(obegin, oit, nbegin, nit); break;
//			default: break;
//			}
//		}
//	}
//
//
////	std::vector<std::vector<esint> > tesize(threads), tenodes(threads);
////	std::vector<std::vector<int> > tetype(threads);
////
////	_fdist.clear();
////	_fdist.reserve(tdistribution.back() + 1);
////	_fdist.push_back(0);
////	for (size_t f = 0; f < fsize.size(); f++) {
////		_fdist.push_back(_fdist.back() + fsize[f]);
////	}
////
////	auto getThreadBegin = [] (const std::vector<esint> &data, const std::vector<esint> &perm, esint eindex) {
////		return std::lower_bound(perm.begin(), perm.end(), eindex, [&] (esint i, esint eindex) {
////			return data[i] < eindex;
////		}) - perm.begin();
////	};
////
////	auto addFaces = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t &triangles, size_t &squares, size_t &index, esint element) {
////		while (index < perm.size() && data[perm[index]] == element) {
////			switch (fsize[perm[index++]]) {
////			case 3: ++triangles; break;
////			case 4:   ++squares; break;
////			}
////		}
////	};
////
////	auto getFace = [&] (const std::vector<esint> &data, const std::vector<esint> &perm, size_t index, esint element, esint n1, esint n2) {
////		while (index < perm.size() && data[perm[index]] == element) {
////			for (esint f = 0; f < fsize[perm[index]]; f++) {
////				if (n1 == fnodes[_fdist[perm[index]] + f] && n2 == fnodes[_fdist[perm[index]] + (f + 1) % fsize[perm[index]]]) {
////					return std::pair<size_t, esint>(index, f);
////				}
////			}
////			++index;
////		}
////		return std::pair<size_t, esint>(index, -1);
////	};
////
////	auto getUnknown = [&] (esint *kbegin, esint *kend, esint *ubegin, esint *uend) {
////		for (auto i = ubegin, j = kbegin; i != uend; ++i) {
////			for (j = kbegin; j != kend; ++j) {
////				if (*i == *j) {
////					break;
////				}
////			}
////			if (j == kend) {
////				return *i;
////			}
////		}
////		return (esint)-1;
////	};
////
////	auto findElementWithSize = [&] (const std::vector<esint> &perm, size_t &index, size_t &max, int size) {
////		while (index < max) { // there is at least one owner
////			if (fsize[perm[index]] == size) {
////				break;
////			} else {
////				++index;
////			}
////		}
////	};
////
////	size_t eoffset = _edist[info::mpi::rank];
////
////	#pragma omp parallel for
////	for (size_t t = 0; t < threads; t++) {
////		std::vector<esint> tsize, tnodes;
////		std::vector<int> ttype;
////
////		size_t oindex = getThreadBegin(owner, owner, tdistribution[t] + eoffset);
////		size_t nindex = getThreadBegin(neighbor, neighbor, tdistribution[t] + eoffset);
////
////		for (size_t e = tdistribution[t] + eoffset; e < tdistribution[t + 1] + eoffset; e++) {
////			size_t obegin = oindex, nbegin = nindex;
////			std::pair<size_t, esint> index;
////			size_t triangles = 0, squares = 0;
////			addFaces(owner, owner, triangles, squares, oindex, e);
////			addFaces(neighbor, neighbor, triangles, squares, nindex, e);
////
////			if (squares == 6 && triangles == 0) {
////				size_t ebegin = tnodes.size();
////				ttype.push_back((int)Element::CODE::HEXA8);
////				tsize.push_back(8);
////				tnodes.insert(tnodes.end(), 8, -1);
////				if (obegin < oindex) { // there is at least one owner
////					tnodes[ebegin + 0] = fnodes[_fdist[owner[obegin]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[owner[obegin]] + 1];
////					tnodes[ebegin + 4] = fnodes[_fdist[owner[obegin]] + 3];
////					tnodes[ebegin + 5] = fnodes[_fdist[owner[obegin]] + 2];
////					++obegin;
////				} else {
////					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nbegin]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nbegin]] + 3];
////					tnodes[ebegin + 4] = fnodes[_fdist[neighbor[nbegin]] + 2];
////					tnodes[ebegin + 5] = fnodes[_fdist[neighbor[nbegin]] + 1];
////					++nbegin;
////				}
////
////				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
////				if (index.first < oindex) {
////					tnodes[ebegin + 2] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
////					tnodes[ebegin + 3] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
////				} else {
////					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
////					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
////					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
////				}
////				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 4], tnodes[ebegin + 5]);
////				if (index.first < oindex) {
////					tnodes[ebegin + 6] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
////					tnodes[ebegin + 7] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
////				} else {
////					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin + 5], tnodes[ebegin + 4]);
////					tnodes[ebegin + 6] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
////					tnodes[ebegin + 7] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
////				}
////				continue;
////			}
////
////			if (squares == 0 && triangles == 4) {
////				size_t ebegin = tnodes.size();
////				ttype.push_back((int)Element::CODE::TETRA4);
////				tsize.push_back(4);
////				tnodes.insert(tnodes.end(), 4, -1);
////				if (obegin < oindex) {
////					tnodes[ebegin + 0] = fnodes[_fdist[owner[obegin]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[owner[obegin]] + 2];
////					tnodes[ebegin + 2] = fnodes[_fdist[owner[obegin]] + 1];
////					++obegin;
////				} else {
////					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nbegin]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nbegin]] + 1];
////					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[nbegin]] + 2];
////					++nbegin;
////				}
////
////				if (obegin < oindex) {
////					tnodes[ebegin + 3] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
////							fnodes.data() + _fdist[owner[obegin]], fnodes.data() + _fdist[owner[obegin] + 1]);
////				} else {
////					tnodes[ebegin + 3] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 3,
////							fnodes.data() + _fdist[neighbor[nbegin]], fnodes.data() + _fdist[neighbor[nbegin] + 1]);
////				}
////				continue;
////			}
////
////			if (squares == 3 && triangles == 2) {
////				size_t ebegin = tnodes.size();
////				ttype.push_back((int)Element::CODE::PRISMA6);
////				tsize.push_back(6);
////				tnodes.insert(tnodes.end(), 6, -1);
////				size_t otria = obegin, ntria = nbegin;
////				findElementWithSize(owner, otria, oindex, 3);
////				findElementWithSize(neighbor, ntria, nindex, 3);
////
////				if (otria < oindex) {
////					tnodes[ebegin + 0] = fnodes[_fdist[owner[otria]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[owner[otria]] + 1];
////					tnodes[ebegin + 2] = fnodes[_fdist[owner[otria]] + 2];
////					findElementWithSize(owner, ++otria, oindex, 3);
////				} else {
////					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[ntria]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[ntria]] + 2];
////					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[ntria]] + 1];
////					findElementWithSize(neighbor, ++ntria, nindex, 3);
////				}
////
////				index = getFace(owner, owner, obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
////				if (fsize[owner[index.first]] == 3) {
////					index = getFace(owner, owner, ++obegin, e, tnodes[ebegin + 1], tnodes[ebegin]);
////				}
////				if (index.first < oindex) {
////					tnodes[ebegin + 3] = fnodes[_fdist[owner[index.first]] + (index.second + 2) % fsize[owner[index.first]]];
////					tnodes[ebegin + 4] = fnodes[_fdist[owner[index.first]] + (index.second + 3) % fsize[owner[index.first]]];
////				} else {
////					index = getFace(neighbor, neighbor, nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
////					if (fsize[neighbor[index.first]] == 3) {
////						index = getFace(neighbor, neighbor, ++nbegin, e, tnodes[ebegin], tnodes[ebegin + 1]);
////					}
////					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[index.first]] + (index.second + 3) % fsize[neighbor[index.first]]];
////					tnodes[ebegin + 4] = fnodes[_fdist[neighbor[index.first]] + (index.second + 2) % fsize[neighbor[index.first]]];
////				}
////
////				if (otria < oindex) {
////					tnodes[ebegin + 5] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
////							fnodes.data() + _fdist[owner[otria]], fnodes.data() + _fdist[owner[otria] + 1]);
////				} else {
////					tnodes[ebegin + 5] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 5,
////							fnodes.data() + _fdist[neighbor[ntria]], fnodes.data() + _fdist[neighbor[ntria] + 1]);
////				}
////				continue;
////			}
////
////			if (squares == 1 && triangles == 4) {
////				size_t ebegin = tnodes.size();
////				ttype.push_back((int)Element::CODE::PYRAMID5);
////				tsize.push_back(5);
////				tnodes.insert(tnodes.end(), 5, -1);
////				size_t osquare = obegin, nsquare = nbegin, otria = obegin, ntria = nbegin;
////				findElementWithSize(owner, osquare, oindex, 4);
////				findElementWithSize(neighbor, nsquare, nindex, 4);
////				findElementWithSize(owner, otria, oindex, 3);
////				findElementWithSize(neighbor, ntria, nindex, 3);
////
////				if (osquare < oindex) {
////					tnodes[ebegin + 0] = fnodes[_fdist[owner[osquare]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[owner[osquare]] + 1];
////					tnodes[ebegin + 2] = fnodes[_fdist[owner[osquare]] + 2];
////					tnodes[ebegin + 3] = fnodes[_fdist[owner[osquare]] + 3];
////				} else {
////					tnodes[ebegin + 0] = fnodes[_fdist[neighbor[nsquare]] + 0];
////					tnodes[ebegin + 1] = fnodes[_fdist[neighbor[nsquare]] + 3];
////					tnodes[ebegin + 2] = fnodes[_fdist[neighbor[nsquare]] + 2];
////					tnodes[ebegin + 3] = fnodes[_fdist[neighbor[nsquare]] + 1];
////				}
////
////				if (otria < oindex) {
////					tnodes[ebegin + 4] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
////							fnodes.data() + _fdist[owner[otria]], fnodes.data() + _fdist[owner[otria] + 1]);
////				} else {
////					tnodes[ebegin + 4] = getUnknown(
////							tnodes.data() + ebegin, tnodes.data() + ebegin + 4,
////							fnodes.data() + _fdist[neighbor[ntria]], fnodes.data() + _fdist[neighbor[ntria] + 1]);
////				}
////				continue;
////			}
////			eslog::error("OpenFOAM parser: an unknown element type with '%ld' triangles and '%ld' squares [ID='%ld'].\n", triangles, squares, e);
////		}
////
////		tesize[t].swap(tsize);
////		tenodes[t].swap(tnodes);
////		tetype[t].swap(ttype);
////	}
////
////	for (size_t t = 0; t < threads; t++) {
////		esize.insert(esize.end(), tesize[t].begin(), tesize[t].end());
////		enodes.insert(enodes.end(), tenodes[t].begin(), tenodes[t].end());
////		etype.insert(etype.end(), tetype[t].begin(), tetype[t].end());
////	}
//}
