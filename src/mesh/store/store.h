
#ifndef SRC_MESH_STORE_STORE_H_
#define SRC_MESH_STORE_STORE_H_

#include <cstddef>
#include <string>
#include <vector>
#include <fstream>

#include "esinfo/mpiinfo.h"
#include "basis/containers/point.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"
#include "mesh/element.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct Store {

	template <typename TBoundaries>
	static void storedata(std::ofstream &os, const std::string &head, const serializededata<TBoundaries, Element*> *data)
	{
		if (data == NULL) {
			return;
		}

		os << head << "\n";
		for (auto elem = data->begin(); elem != data->end(); ++elem) {
			os << "[ ";
			for (auto i = elem->begin(); i != elem->end(); ++i) {
				os << (int)(*i)->code << " ";
			}
			os << "]\n";
		}
		os << "\n";
	}

	template <typename TBoundaries>
	static void storedata(std::ofstream &os, const std::string &head, const serializededata<TBoundaries, Point> *data)
	{
		if (data == NULL) {
			return;
		}

		os << head << "\n";
		for (auto elem = data->begin(); elem != data->end(); ++elem) {
			os << "[ ";
			for (auto i = elem->begin(); i != elem->end(); ++i) {
				os << "<" << i->x << " " << i->y << " " << i->z << "> ";
			}
			os << "]\n";
		}
		os << "\n";
	}

	template <typename TBoundaries, typename TData>
	static void storedata(std::ofstream &os, const std::string &head, const serializededata<TBoundaries, TData> *data)
	{
		if (data == NULL) {
			return;
		}

		os << head << "\n";
		for (auto elem = data->begin(); elem != data->end(); ++elem) {
			os << "[ ";
			for (auto i = elem->begin(); i != elem->end(); ++i) {
				os << *i << " ";
			}
			os << "]\n";
		}
		os << "\n";
	}

	static std::vector<esint> gatherDistribution(esint size)
	{
		std::vector<esint> result(info::mpi::size + 1);
		Communication::allGather(&size, result.data(), MPITools::getType(size).mpisize, MPITools::getType(size).mpitype);
		utils::sizesToOffsets(result);
		return result;
	}

	template <typename TType>
	static std::vector<TType> gatherDistribution(std::vector<TType> &distribution, TType offset)
	{
		std::vector<TType> odistribution = distribution;
		for (size_t i = 0; i < odistribution.size(); i++) {
			odistribution[i] += offset;
		}
		std::vector<TType> result;
		Communication::gatherUnknownSize(odistribution, result);
		utils::removeDuplicates(result);
		Communication::broadcastUnknownSize(result);

		return result;
	}
};

}


#endif /* SRC_MESH_STORE_STORE_H_ */
