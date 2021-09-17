
#ifndef SRC_ANALYSIS_SCHEME_SCHEME_H_
#define SRC_ANALYSIS_SCHEME_SCHEME_H_

namespace espreso {

template <typename T>
static void clear(T *t)
{
	if (t) {
		delete t;
	}
}

template <typename T, typename ...Other>
static void clear(T *t, Other... other)
{
	clear(t);
	clear(other...);
}

}

#endif /* SRC_ANALYSIS_SCHEME_SCHEME_H_ */
