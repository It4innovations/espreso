
#ifndef SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_
#define SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_

#include "edata.h"
#include "tarray.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/packing.h"

#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData>
class serializededata {

public:
	static void balance(size_t esize, std::vector<std::vector<TEData> > &data, const std::vector<size_t> *distribution = NULL)
	{
//		profiler::syncstart("balance_fixed");
		std::vector<size_t> _distribution;
		if (distribution == NULL) {
			size_t size = 0;
			for (size_t t = 0; t < data.size(); t++) {
				size += data[t].size();
			}

			_distribution = tarray<size_t>::distribute(data.size(), size);
		} else {
			_distribution = *distribution;
		}

//		profiler::syncparam("size", _distribution.back());
		for (size_t t = 0, tt = 0; t < data.size(); tt = ++t) {
			while (++tt && data[t].size() < esize * (_distribution[t + 1] - _distribution[t])) {
				size_t diff = esize * (_distribution[t + 1] - _distribution[t]) - data[t].size();
				if (diff < data[tt].size()) {
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].begin() + diff);
					data[tt].erase(data[tt].begin(), data[tt].begin() + diff);
				} else {
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].end());
					data[tt].clear();

				}
			}
			if (data[t].size() > esize * (_distribution[t + 1] - _distribution[t])) {
				size_t diff = data[t].size() - esize * (_distribution[t + 1] - _distribution[t]);
				data[t + 1].insert(data[t + 1].begin(), data[t].end() - diff, data[t].end());
				data[t].erase(data[t].end() - diff, data[t].end());
			}
		}
//		profiler::syncend("balance_fixed");
	}

	static void balance(std::vector<std::vector<TEBoundaries> > &boundaries, std::vector<std::vector<TEData> > &data, const std::vector<size_t> *distribution = NULL)
	{
//		profiler::syncstart("balance_various");
		size_t size = 0;
		std::vector<size_t> sizes(boundaries.size());
		for (size_t t = 0; t < boundaries.size(); t++) {
			sizes[t] = boundaries[t].size();
			size += boundaries[t].size();
		}
		--sizes[0];
//		profiler::syncparam("size", size);
		if (size == 0) {
//			profiler::syncend("balance_various");
			return;
		}

		std::vector<size_t> _distribution;
		if (distribution == NULL) {
			_distribution = tarray<size_t>::distribute(data.size(), size - 1);
		} else {
			_distribution = *distribution;
		}

		for (size_t t = 0, tt = 0; t < boundaries.size(); tt = ++t) {
			while (++tt && sizes[t] < _distribution[t + 1] - _distribution[t]) {
				size_t diff = _distribution[t + 1] - _distribution[t] - sizes[t];
				if (diff < sizes[tt]) {
					size_t ttt = 0;
					while (boundaries[tt - ++ttt].size() == 0);
					size_t ediff = *(boundaries[tt].begin() + diff - 1) - boundaries[tt - ttt].back();
					data[t].insert(data[t].end(), data[tt].begin(), data[tt].begin() + ediff);
					data[tt].erase(data[tt].begin(), data[tt].begin() + ediff);

					boundaries[t].insert(boundaries[t].end(), boundaries[tt].begin(), boundaries[tt].begin() + diff);
					boundaries[tt].erase(boundaries[tt].begin(), boundaries[tt].begin() + diff);
					sizes[t]  += diff;
					sizes[tt] -= diff;
				} else {
					sizes[t]  += boundaries[tt].size();
					sizes[tt] -= boundaries[tt].size();

					data[t].insert(data[t].end(), data[tt].begin(), data[tt].end());
					data[tt].clear();

					boundaries[t].insert(boundaries[t].end(), boundaries[tt].begin(), boundaries[tt].end());
					boundaries[tt].clear();
				}
			}
			if (sizes[t] > _distribution[t + 1] - _distribution[t]) {
				size_t diff = sizes[t] - (_distribution[t + 1] - _distribution[t]);
				size_t ediff = boundaries[t].back() - *(boundaries[t].end() - diff - 1);
				data[t + 1].insert(data[t + 1].begin(), data[t].end() - ediff, data[t].end());
				data[t].erase(data[t].end() - ediff, data[t].end());

				boundaries[t + 1].insert(boundaries[t + 1].begin(), boundaries[t].end() - diff, boundaries[t].end());
				boundaries[t].erase(boundaries[t].end() - diff, boundaries[t].end());
				sizes[t]     -= diff;
				sizes[t + 1] += diff;
			}
		}
//		profiler::syncend("balance_various");
	}

private:
	template<class TIterator, typename TIteratorEData>
	class iterator_base {

	public:
//		bool operator< (const TIterator &other) const { return _edata._begin <  other._edata._begin; }
//		bool operator> (const TIterator &other) const { return _edata._begin >  other._edata._begin; }
//		bool operator<=(const TIterator &other) const { return _edata._begin <= other._edata._begin; }
//		bool operator>=(const TIterator &other) const { return _edata._begin >= other._edata._begin; }
		bool operator==(const TIterator &other) const { return _element == other._element && _edata._begin == other._edata._begin; }
		bool operator!=(const TIterator &other) const { return _element != other._element || _edata._begin != other._edata._begin; }

		TIterator& operator++() { return move( 1); }
		TIterator& operator--() { return move(-1); }
		TIterator  operator++(int) {TIterator tmp(*static_cast<TIterator*>(this)); operator++(); return tmp; }
		TIterator  operator--(int) {TIterator tmp(*static_cast<TIterator*>(this)); operator--(); return tmp; }
		template <typename TType> TIterator  operator+ (TType n) { return TIterator(*static_cast<TIterator*>(this)).move( n); }
		template <typename TType> TIterator  operator- (TType n) { return TIterator(*static_cast<TIterator*>(this)).move(-n); }
		template <typename TType> TIterator& operator+=(TType n) { return move( n); }
		template <typename TType> TIterator& operator-=(TType n) { return move(-n); }

	protected:
		iterator_base(const TEBoundaries *begin, const TEBoundaries *element, const TEBoundaries *end, TIteratorEData *edata)
		: _element(begin), _end(end), _edata(edata, edata)
		{
			move(element - begin);
		}

		iterator_base(size_t edatasize, TIteratorEData *edata)
		: _element(NULL), _end(NULL), _edata(edata, edata)
		{
			_edata._end += edatasize;
		}

		template <typename TType>
		TIterator& move(TType n)
		{
			if (_element == NULL) {
				size_t size = _edata._end - _edata._begin;
				_edata._begin += n * size;
				_edata._end += n * size;
			} else {
				_edata._begin += *(_element + n) - *(_element);
				_element += n;
				if (_element != _end) {
					_edata._end = _edata._begin + *(_element + 1) - *_element;
				} else {
					_edata._end = _edata._begin;
				}
			}
			return static_cast<TIterator&>(*this);
		}

		const TEBoundaries* _element;
		const TEBoundaries* _end;
		edata<TIteratorEData> _edata;
	};

public:

	class iterator: public iterator_base<iterator, TEData>
	{
		friend class serializededata<TEBoundaries, TEData>;
	public:
		edata<TEData>& operator*()  { return  this->_edata; }
		edata<TEData>* operator->() { return &this->_edata; }

	private:
		iterator(TEBoundaries *begin, TEBoundaries *element, TEBoundaries *end, TEData *edata)
		: iterator_base<iterator, TEData>(begin, element, end, edata) { }
		iterator(size_t edatasize, TEData *edata)
		: iterator_base<iterator, TEData>(edatasize, edata) { }
	};

	class const_iterator: public iterator_base<const_iterator, const TEData>
	{
		friend class serializededata<TEBoundaries, TEData>;
	public:

		edata<const TEData>& operator*()  { return  this->_edata; }
		edata<const TEData>* operator->() { return &this->_edata; }

	private:
		const_iterator(const TEBoundaries *begin, const TEBoundaries *element, const TEBoundaries *end, const TEData *edata)
		: iterator_base<const_iterator, const TEData>(begin, element, end, edata) { }
		const_iterator(size_t edatasize, TEData *edata)
		: iterator_base<const_iterator, const TEData>(edatasize, edata) { }
	};

	// data are uniform
	serializededata(size_t edatasize, tarray<TEData> &&edata)
	: _eboundaries(0, 0), _edata(std::move(edata)), _edatasize(edatasize) { inititerators(edatasize); }

	serializededata(size_t edatasize, const std::vector<size_t> &distribution, TEData init = TEData{})
	: _eboundaries(0, 0), _edata(distribution, edatasize, init), _edatasize(edatasize) { inititerators(edatasize); }

	// data are non-uniform
	serializededata(tarray<TEBoundaries> &&eboundaries, tarray<TEData> &&edata)
	: _eboundaries(std::move(eboundaries)), _edata(std::move(edata)), _edatasize(-1) { inititerators(); }

	serializededata(const serializededata<TEBoundaries, TEData> &other)
	: _eboundaries(other._eboundaries), _edata(other._edata), _edatasize(other._edatasize) { _edatasize != -1 ? inititerators(_edatasize) : inititerators(); }

	serializededata(serializededata<TEBoundaries, TEData> &&other)
	: _eboundaries(std::move(other._eboundaries)), _edata(std::move(other._edata)), _edatasize(std::move(other._edatasize)) { inititerators(); }

	serializededata<TEBoundaries, TEData>& operator=(const serializededata<TEBoundaries, TEData> &other)
	{
		if (this != &other) {
			_eboundaries = other._eboundaries;
			_edata = other._edata;
			_edatasize = other._edatasize;
			inititerators();
		}
		return *this;
	}
	serializededata<TEBoundaries, TEData>& operator=(serializededata<TEBoundaries, TEData> &&other)
	{
		if (this != &other) {
			_eboundaries = std::move(other._eboundaries);
			_edata = std::move(other._edata);
			_edatasize = std::move(other._edatasize);
			inititerators();
		}
		return *this;
	}

	size_t threads() const { return _edata.threads(); }
	size_t structures() const
	{
		if (_eboundaries.size()) {
			return _eboundaries.size() - 1;
		} else {
			return _edata.size() / _edatasize;
		}
	}

	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
	{
		if (_eboundaries.size()) {
			permuteNonUniformData(permutation, distribution);
		} else {
			permuteUniformData(permutation, distribution);
		}
	}

	iterator       begin()        { return _iterator.front(); }
	const_iterator begin()  const { return _constiterator.front(); }
	const_iterator cbegin() const { return _constiterator.front(); }
	iterator       end()          { return _iterator.back(); }
	const_iterator end()    const { return _constiterator.back(); }
	const_iterator cend()   const { return _constiterator.back(); }

	iterator       begin (size_t thread)       { return _iterator[thread]; }
	const_iterator begin (size_t thread) const { return _constiterator[thread]; }
	const_iterator cbegin(size_t thread) const { return _constiterator[thread]; }
	iterator       end   (size_t thread)       { return _iterator[thread + 1]; }
	const_iterator end   (size_t thread) const { return _constiterator[thread + 1]; }
	const_iterator cend  (size_t thread) const { return _constiterator[thread + 1]; }

	tarray<TEBoundaries>&       boundarytarray()       { return _eboundaries; }
	const tarray<TEBoundaries>& boundarytarray() const { return _eboundaries; }
	tarray<TEData>&             datatarray()           { return _edata; }
	const tarray<TEData>&       datatarray()     const { return _edata; }

	size_t packedSize() const
	{
		return sizeof(esint) + _eboundaries.packedSize() + _edata.packedSize();
	}

	void pack(char* &p) const
	{
		memcpy(p, &_edatasize, sizeof(esint));
		p += sizeof(esint);

		_eboundaries.pack(p);
		_edata.pack(p);
	}

	void unpack(const char* &p)
	{
		memcpy(&_edatasize, p, sizeof(esint));
		p += sizeof(esint);

		_eboundaries.unpack(p);
		_edata.unpack(p);

		if (_edatasize == -1) {
			inititerators();
		} else {
			inititerators(_edatasize);
		}
	}

private:
	void inititerators()
	{
		_iterator = std::vector<iterator>(threads() + 1, iterator(_eboundaries.begin(), _eboundaries.begin(), _eboundaries.end() - 1, _edata.begin()));
		_constiterator = std::vector<const_iterator>(threads() + 1, const_iterator(_eboundaries.begin(), _eboundaries.begin(), _eboundaries.end() - 1, _edata.begin()));

		for (size_t t = 1; t <= threads(); t++) {
			_iterator[t] += _eboundaries.distribution()[t] - 1;
			_constiterator[t] += _eboundaries.distribution()[t] - 1;
		}
	}

	void inititerators(size_t edatasize)
	{
		_iterator = std::vector<iterator>(threads() + 1, iterator(edatasize, _edata.begin()));
		_constiterator = std::vector<const_iterator>(threads() + 1, const_iterator(edatasize, _edata.begin()));

		for (size_t t = 1; t <= threads(); t++) {
			_iterator[t] += _edata.distribution()[t] / edatasize;
			_constiterator[t] += _edata.distribution()[t] / edatasize;
		}
	}

	void permuteUniformData(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
	{
//		profiler::syncstart("permute_uniform_data");
//		profiler::syncparam("esize", _edatasize);
//		profiler::syncparam("elements", permutation.size());
		std::vector<std::vector<TEData> > pdata(distribution.size() - 1);

		#pragma omp parallel for
		for (size_t t = 0; t < distribution.size() - 1; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
				pdata[t].insert(pdata[t].end(), _edata.data() + _edatasize * permutation[i], _edata.data() + _edatasize * (permutation[i] + 1));
			}
		}

		_edata = tarray<TEData>(pdata);
		inititerators(_edatasize);
//		profiler::syncend("permute_uniform_data");
	}

	void permuteNonUniformData(const std::vector<esint> &permutation, const std::vector<size_t> &providedDistribution)
	{
//		profiler::syncstart("permute_non_uniform_data");
//		profiler::syncparam("elements", permutation.size());
//		profiler::syncparam("datasize", _eboundaries.back());
		std::vector<std::vector<TEBoundaries> > pboundaries(providedDistribution.size() - 1);
		std::vector<std::vector<TEData> > pdata(providedDistribution.size() - 1);
		std::vector<size_t> distribution = providedDistribution;
		if (_eboundaries.distribution().back() > distribution.back()) {
			for (size_t t = 1; t < providedDistribution.size(); t++) {
				distribution[t] += 1;
			}
		}

		pboundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < providedDistribution.size() - 1; t++) {
			for (size_t e = (t == 0 ? 1 : distribution[t]); e < distribution[t + 1]; ++e) {
				pboundaries[t].push_back(_eboundaries.data()[permutation[e - 1] + 1] - _eboundaries.data()[permutation[e - 1]]);
				pdata[t].insert(
						pdata[t].end(),
						_edata.data() + _eboundaries.data()[permutation[e - 1]],
						_edata.data() + _eboundaries.data()[permutation[e - 1] + 1]);
				if (pboundaries[t].size() > 1) {
					pboundaries[t].back() += *(pboundaries[t].end() - 2);
				}
			}
		}

		std::vector<size_t> offsets;
		for (size_t t = 0; t < pboundaries.size(); t++) {
			offsets.push_back(pboundaries[t].size() ? pboundaries[t].back() : 0);
		}

		TEBoundaries sum = 0;
		for (size_t i = 0; i < offsets.size(); i++) {
			TEBoundaries tmp = offsets[i];
			offsets[i] = sum;
			sum += tmp;
		}

		#pragma omp parallel for
		for (size_t t = 0; t < pboundaries.size(); t++) {
			size_t offset = offsets[t];
			for (size_t i = 0; i < pboundaries[t].size(); i++) {
				pboundaries[t][i] += offset;
			}
		}

		_eboundaries = tarray<TEBoundaries>(pboundaries);
		_edata = tarray<TEData>(pdata);
		inititerators();
//		profiler::syncend("permute_non_uniform_data");
	}

	tarray<TEBoundaries> _eboundaries;
	tarray<TEData> _edata;

	esint _edatasize;
	std::vector<iterator> _iterator;
	std::vector<const_iterator> _constiterator;
};

namespace utils {

template <typename TEBoundaries, typename TEData>
inline size_t packedSize(serializededata<TEBoundaries, TEData> *data)
{
	if (data != NULL) {
		return data->packedSize() + 1 + 3 * sizeof(size_t);
	}
	return 1;
}

template <typename TEBoundaries, typename TEData>
inline void pack(serializededata<TEBoundaries, TEData> *data, char* &p)
{
	pack(data != NULL, p);
	if (data != NULL) {
		pack(data->threads(), p);
		pack(data->boundarytarray().size(), p);
		pack(data->datatarray().size(), p);
		data->pack(p);
	}
}

template <typename TEBoundaries, typename TEData>
inline void unpack(serializededata<TEBoundaries, TEData> *&data, const char* &p)
{
	if (data != NULL) {
		delete data;
		data = NULL;
	}

	bool notnull;
	unpack(notnull, p);
	if (notnull) {
		size_t threads, bsize, dsize;
		unpack(threads, p);
		unpack(bsize, p);
		unpack(dsize, p);
		data = new serializededata<TEBoundaries, TEData>(tarray<TEBoundaries>(threads, bsize), tarray<TEData>(threads, dsize));
		data->unpack(p);
	}
}

} // namespace utils
} // namespace espreso


#endif /* SRC_BASIS_CONTAINERS_SERIALIZEDEDATA_H_ */
