
#ifndef SRC_ESINFO_STEPINFO_H_
#define SRC_ESINFO_STEPINFO_H_

namespace espreso {
namespace step {

	enum class TYPE: int {
		TIME = 1,
		FREQUENCY = 2,
		FTT = 4
	};

	inline TYPE operator|(TYPE t1, TYPE t2) { return static_cast<TYPE>(static_cast<int>(t1) | static_cast<int>(t2)); }
	inline TYPE operator&(TYPE t1, TYPE t2) { return static_cast<TYPE>(static_cast<int>(t1) & static_cast<int>(t2)); }

	extern int loadstep;
	extern int substep;
	extern int iteration;

	extern TYPE type;

namespace duplicate {
	extern int instances;
	extern int offset;
	extern int size;
	extern int totalsize;
}

namespace time {
	extern double start;
	extern double current;
	extern double shift; // difference between current and previous time
	extern double final;
	extern double precision;

	inline bool isLast() {
		if (duplicate::instances == 1) {
			return current == final;
		} else {
			return substep - duplicate::offset == duplicate::size;
		}
	}
}

namespace frequency {
	extern double start;
	extern double current;
	extern double shift; // difference between current and previous frequency
	extern double final;
	extern double precision;
	extern double angular;

	inline bool isLast() {
		if (duplicate::instances == 1) {
			return current == final;
		} else {
			return substep - duplicate::offset == duplicate::size;
		}
	}
}

namespace ftt { // frequency to time
	extern int step;
	extern int steps;

	extern double time;
	extern double period;

	inline bool isFirst() { return step == 0; }
	inline bool isLast() { return step + 1 == steps; }
}

inline bool isInitial()
{
	return loadstep == 0 && substep == 0 && iteration == 0;
}

inline bool isFirst()
{
	switch (type) {
	case TYPE::TIME:
	case TYPE::FREQUENCY:
		return substep == 0 && iteration == 0;
	case TYPE::FTT:
		return ftt::isFirst();
	}
	return false;
}

inline bool isLast()
{
	switch (type) {
	case TYPE::TIME:
		return time::isLast();
	case TYPE::FREQUENCY:
		return frequency::isLast();
	case TYPE::FTT:
		return ftt::isLast();
	}
	return false;
}

};

}


#endif /* SRC_ESINFO_STEPINFO_H_ */
