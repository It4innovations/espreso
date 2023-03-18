
#include "tableinterpolationevaluator.h"

#include "esinfo/eslog.hpp"

using namespace espreso;

TableInterpolationEvaluator::TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table)
: _table(table)
{
	if (!table.size()) {
		eslog::globalerror("Interpolation table with zero size.\n");
	}
}

double TableInterpolationEvaluator::evaluate(int t) const
{
//	for (esint i = 0; i < size; ++i) {
//		if (params._temp[i] < _table[0].first) {
//			results[i * increment] = _table[0].second;
//			break;
//		}
//		for (size_t j = 0; j < _table.size() - 1; j++) {
//			if (_table[j].first < params._temp[i] && params._temp[i] < _table[j + 1].first) {
//				double a = _table[j].first, b = _table[j + 1].first;
//				double va = _table[j].second, vb = _table[j + 1].second;
//				results[i * increment] = va + (vb - va) * (params._temp[i] - a) / (b - a);
//				break;
//			}
//		}
//		results[i * increment] = _table.back().second;
//	}
	return 0;
}

