
#ifndef SRC_INPUT_PARSERS_MIXEDELEMENTSPARSER_H_
#define SRC_INPUT_PARSERS_MIXEDELEMENTSPARSER_H_

#include <cstddef>
#include <vector>
#include <functional>

namespace espreso {

struct InputFile;

class MixedElementsParser {
public:
    struct ElementData {
        esint *data;
        size_t size;
    };

    void add(esint *data, size_t size);
    void parse(std::function<esint(size_t, esint)> enodes);

    std::vector<int> invalid;
    std::vector<esint> first, missing, offset, nelements;
protected:
    void recognize(size_t index, esint &esize, esint &coffset, esint &missing, std::function<esint(size_t, esint)> enodes);
    void fix(std::vector<esint> &lmissing, std::vector<esint> &esize, std::vector<int> &last, std::function<esint(size_t, esint)> enodes);

    std::vector<ElementData> _elements;
};

}

#endif /* SRC_INPUT_PARSERS_MIXEDELEMENTSPARSER_H_ */
