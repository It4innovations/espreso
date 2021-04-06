
#ifndef SRC_ESINFO_SYSTEMINFO_H_
#define SRC_ESINFO_SYSTEMINFO_H_

namespace espreso {
namespace info {
namespace system {

	enum class OPERATIONSYSTEM {
		UNIX
	};

	enum class INSTRUCTIONSET {
		SSE,
		AVX
	};

	enum class BUILD_MODE {
		RELEASE,
		DEVEL,
		DEBUG,
		PROFILE
	};

	// dummy yet
	OPERATIONSYSTEM os();
	INSTRUCTIONSET instructionSet();

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
	int MPI_Init_thread_provided_level();

	void setSignals();
}
}
}


#endif /* SRC_ESINFO_SYSTEMINFO_H_ */
