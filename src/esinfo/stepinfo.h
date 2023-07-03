
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

struct Step {
	int loadstep = 0, loadsteps = 1;
	int substep = 0, substeps = 1;
	int iteration = 0, iterations = 1;

	TYPE type = TYPE::TIME;
} extern step, outstep;

struct Duplicate {
	int instances = 1;
	int offset = 0;
	int size = 1;
	int totalsize = 1;
} extern duplicate, outduplicate;

struct Time {
	double start = 0;
	double current = 0;
	double shift = 1; // difference between current and previous time
	double final = 1;
	double precision = 1e-8;

	inline bool isLast(const Step &step, const Duplicate &duplicate) const {
		if (duplicate.instances == 1) {
			return current == final;
		} else {
			return step.substep - duplicate.offset == duplicate.size;
		}
	}
} extern time, outtime;

struct Frequency {
	double start = 0;
	double current = 0;
	double shift = 1; // difference between current and previous frequency
	double final = 1;
	double precision = 1e-8;
	double angular = 0;

	inline bool isLast(const Step &step, const Duplicate &duplicate) const {
		if (duplicate.instances == 1) {
			return current == final;
		} else {
			return step.substep - duplicate.offset == duplicate.size;
		}
	}
} extern frequency, outfrequency;

struct Ftt { // frequency to time
	int step = 0;
	int steps = 1;

	double time = 0;
	double period = 1;

	inline bool isFirst() const { return step == 0; }
	inline bool isLast() const { return step + 1 == steps; }
} extern ftt, outftt;


inline bool isInitial(const Step &step = step)
{
	return step.loadstep == 0 && step.substep == 0 && step.iteration == 0;
}

inline bool isFirst(const Step &step = step, const Ftt &ftt = ftt)
{
	switch (step.type) {
	case TYPE::TIME:
	case TYPE::FREQUENCY:
		return step.substep == 0 && step.iteration == 0;
	case TYPE::FTT:
		return ftt.isFirst();
	}
	return false;
}

inline bool isLast(const Step &step = step, const Duplicate &duplicate = duplicate, const Time &time = time, const Frequency &frequency = frequency, const Ftt &ftt = ftt)
{
	switch (step.type) {
	case TYPE::TIME:
		return time.isLast(step, duplicate);
	case TYPE::FREQUENCY:
		return frequency.isLast(step, duplicate);
	case TYPE::FTT:
		return ftt.isLast();
	}
	return false;
}

inline void toOut()
{
	outstep = step;
	outduplicate = duplicate;
	outtime = time;
	outfrequency = frequency;
	outftt = ftt;
}

};

}


#endif /* SRC_ESINFO_STEPINFO_H_ */
