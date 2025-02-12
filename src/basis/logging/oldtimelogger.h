
#ifndef SRC_BASIS_LOGGING_OLDTIMELOGGER_H_
#define SRC_BASIS_LOGGING_OLDTIMELOGGER_H_

// TODO: dummy only for assign verebosity and m parameter

#include "verbosity.h"

namespace espreso {

// only for verbosity parameter -m
class OldTimeLogger: public Verbosity<OldTimeLogger, 'm'> {

public:
    void initOutput()
    {

    }

    void start(const char* region, const char* section)
    {

    }

    void checkpoint(const char* region)
    {

    }

    void accumulated(const char* region)
    {

    }

    void end(const char* region)
    {

    }

    void param(const char* name, const int &value)
    {

    }

    void param(const char* name, const long &value)
    {

    }

    void param(const char* name, const long unsigned int &value)
    {

    }

    void param(const char* name, const double &value)
    {

    }

    void param(const char* name, const char* value)
    {

    }

    void ln()
    {

    }

    void nextLoadStep(int step)
    {

    }

    void output(const char* msg, VerboseArg::COLOR color)
    {

    }

    void error(const char* msg)
    {

    }

    void finish()
    {

    }
};

}



#endif /* SRC_BASIS_LOGGING_OLDTIMELOGGER_H_ */
