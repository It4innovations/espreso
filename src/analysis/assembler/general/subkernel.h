
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_SUBKERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_SUBKERNEL_H_


namespace espreso {

struct SubKernel {
    enum Action: int {
        VOID       = 1 << 0,
        PREPROCESS = 1 << 1,
        ASSEMBLE   = 1 << 2,
        REASSEMBLE = 1 << 3,
        ITERATION  = 1 << 4,
        SOLUTION   = 1 << 5
    };

    bool isconst, isactive;
    bool needCoordinates, needTemperature;
    Action action;

    SubKernel(): isconst(1), isactive(0), needCoordinates(false), needTemperature(false), action(VOID) {}

    void setActiveness(Action action)
    {
        setActiveness(action, true);
    }

    void setActiveness(Action action, int guard)
    {
        isactive = isactive && (this->action & action) && guard;
    }
};

inline SubKernel::Action  operator| (SubKernel::Action  a1, SubKernel::Action a2) { return static_cast<SubKernel::Action>(static_cast<int>(a1) | static_cast<int>(a2)); }
inline SubKernel::Action  operator& (SubKernel::Action  a1, SubKernel::Action a2) { return static_cast<SubKernel::Action>(static_cast<int>(a1) & static_cast<int>(a2)); }
inline SubKernel::Action& operator|=(SubKernel::Action &a1, SubKernel::Action a2) { a1 = a1 | a2; return a1; }
inline SubKernel::Action& operator&=(SubKernel::Action &a1, SubKernel::Action a2) { a1 = a1 & a2; return a1; }

}


#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_SUBKERNEL_H_ */
