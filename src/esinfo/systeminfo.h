
#ifndef SRC_ESINFO_SYSTEMINFO_H_
#define SRC_ESINFO_SYSTEMINFO_H_

namespace espreso {
namespace info {
namespace system {

	enum class OPERATIONSYSTEM {
		UNIX
	};

	enum class BUILD_MODE {
		RELEASE,
		DEVEL,
		DEBUG,
		PROFILE
	};

	struct CPUInfo {
		char modelName[200];
		int sockets = 0;
		int cores;
		bool hyperthreading = 0;
	};

	// dummy yet
	OPERATIONSYSTEM os();

	constexpr BUILD_MODE buildmode()
	{
		#ifdef __ESMODE__
			return BUILD_MODE::__ESMODE__;
		#else
			return BUILD_MODE::RELEASE;
		#endif
	}

	const char* commit();
	const char* cxx();
	const char* buildpath();
	const char* cxxflags();
	const char* simd();

	int MPI_Init_thread_provided_level();

	CPUInfo cpuinfo();
	int hwthreads();
	long mpiAffinity();
	void forceAffinity(int size);
	bool pinningIntersection();
	long memoryTotal();
	long memoryAvail();

	void setSignals();
	void print();
}
}
}


#endif /* SRC_ESINFO_SYSTEMINFO_H_ */
