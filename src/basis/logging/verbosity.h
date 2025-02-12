
#ifndef SRC_BASIS_LOGGING_VERBOSITY_H_
#define SRC_BASIS_LOGGING_VERBOSITY_H_

namespace espreso {

struct VerboseArg {
    enum class COLOR {
        WHITE,
        RED,
        GREEN,
        YELLOW,
        BLUE,
        MAGENTA,
        CYAN
    };

    char argflag;
    int level;
    int verbosity;
    int finishing;
    int always;
    int rank;
    int grank;

    VerboseArg(char argflag)
    : argflag(argflag),
      level(0), verbosity(1),
      finishing(0), always(0),
      rank(0), grank(0) {}

    bool isAllowed() const {
        return always || level <= verbosity;
    }
};

template <typename T, char C>
struct Verbosity: public VerboseArg {
    Verbosity(): VerboseArg(C) {}
};

}



#endif /* SRC_BASIS_LOGGING_VERBOSITY_H_ */
