#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <cstdlib>
#include <cstdio>

class Settings
{

public:

	enum FLAGS {
		DYNAMIC,
		FLAGS_SIZE
	};

	static void setFlag(int flag, bool value)
	{
		if (flag < 0 || flag >= FLAGS_SIZE) {
			fprintf(stderr, "Unknown flag.\n");
			exit(EXIT_FAILURE);
		}
		_flags[flag] = value;
	}

	static bool getFlag(int flag)
	{
		if (flag < 0 || flag >= FLAGS_SIZE) {
			fprintf(stderr, "Unknown flag.\n");
			exit(EXIT_FAILURE);
		}
		return _flags[flag];
	}


private:
	static bool _flags[FLAGS_SIZE];
};



#endif /* SETTINGS_H_ */
