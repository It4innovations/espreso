
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
};

struct Duplicate {
	int instances = 1;
	int offset = 0;
	int size = 1;
	int totalsize = 1;
};

struct Time {
	double start = 0;
	double current = 0;
	double shift = 1; // difference between current and previous time
	double final = 1;
	double precision = 1e-8;
};

struct Frequency {
	double start = 0;
	double current = 0;
	double shift = 1; // difference between current and previous frequency
	double final = 1;
	double precision = 1e-8;
	double angular = 0;
};

struct Ftt { // frequency to time
	int step = 0;
	int steps = 1;

	double time = 0;
	double period = 1;
};

inline bool isLast(const Step &step)
{
	return step.substep + 1 == step.substeps;
}

};

}


#endif /* SRC_ESINFO_STEPINFO_H_ */
