
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

	int isconst, update;
	Action action;

	ActionOperator(): isconst(1), update(1), action(Action::VOID) {}
	virtual ~ActionOperator() {}

	virtual void move(int n) {};
};

inline ActionOperator::Action  operator| (ActionOperator::Action  a1, ActionOperator::Action a2) { return static_cast<ActionOperator::Action>(static_cast<int>(a1) | static_cast<int>(a2)); }
inline ActionOperator::Action  operator& (ActionOperator::Action  a1, ActionOperator::Action a2) { return static_cast<ActionOperator::Action>(static_cast<int>(a1) & static_cast<int>(a2)); }
inline ActionOperator::Action& operator|=(ActionOperator::Action &a1, ActionOperator::Action a2) { a1 = a1 | a2; return a1; }
inline ActionOperator::Action& operator&=(ActionOperator::Action &a1, ActionOperator::Action a2) { a1 = a1 & a2; return a1; }

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
