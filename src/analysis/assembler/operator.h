
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

namespace espreso {

struct ActionOperator {
	enum Action: int {
		VOID       = 1 << 0,
		ASSEMBLE   = 1 << 1,
		REASSEMBLE = 1 << 2,
		FILL       = 1 << 3,
		SOLUTION   = 1 << 4
	};

	static inline void removeSolution(Action &action);

	int isconst, update;
	Action action;

	ActionOperator(): isconst(1), update(1), action(Action::VOID) {}
	virtual ~ActionOperator() {}

	virtual void setTime(double time, int t) {};
	virtual void setFrequency(double frequency, int t) {};
	virtual void move(int n) {};

	virtual const char* name() const =0;
};

inline ActionOperator::Action  operator| (ActionOperator::Action  a1, ActionOperator::Action a2) { return static_cast<ActionOperator::Action>(static_cast<int>(a1) | static_cast<int>(a2)); }
inline ActionOperator::Action  operator& (ActionOperator::Action  a1, ActionOperator::Action a2) { return static_cast<ActionOperator::Action>(static_cast<int>(a1) & static_cast<int>(a2)); }
inline ActionOperator::Action& operator|=(ActionOperator::Action &a1, ActionOperator::Action a2) { a1 = a1 | a2; return a1; }
inline ActionOperator::Action& operator&=(ActionOperator::Action &a1, ActionOperator::Action a2) { a1 = a1 & a2; return a1; }

inline void ActionOperator::removeSolution(Action &action) { action &= Action::ASSEMBLE | Action::REASSEMBLE | Action::FILL; }

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
