
#ifndef SRC_BASIS_UTILITIES_COMMUNICATION_H_
#define SRC_BASIS_UTILITIES_COMMUNICATION_H_

namespace espreso {

struct Communication {

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer);

	template <typename Ttype>
	static bool gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets);
};


}

#include "communication.hpp"




#endif /* SRC_BASIS_UTILITIES_COMMUNICATION_H_ */
