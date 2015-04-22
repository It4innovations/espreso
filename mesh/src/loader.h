#ifndef LOADING_H_
#define LOADING_H_

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

class Loader
{

public:
	static size_t getLinesCount(const char* fileName);

private:

	struct TestEOL {
		bool operator()(char c)
		{
			return c == '\n';
		}
	};
};


#endif /* LOADING_H_ */
