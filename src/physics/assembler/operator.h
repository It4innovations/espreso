
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

namespace espreso {

#define GET_NAME(structure) const char* name() { return #structure; }

class HeatTransferKernelOpt;

struct Operator {
	static const int print = 0;

	class Link {
		int interval, version, update, isset;

	public:
		Link(int interval): interval(interval), version(0), update(0), isset(0) {}

		operator bool() const { return update; }

		Link& inputs(int version)
		{
			this->version = version;
			isset = 1;
			return *this;
		}

		template <class TParameter>
		Link& inputs(const TParameter &data)
		{
			if (version < data.version[interval]) {
				version = data.version[interval];
			}
			if (data.isset) {
				isset = data.isset;
			}
			return *this;
		}

		template <class TParameter, class ...Other>
		Link& inputs(const TParameter &data, const Other&...other)
		{
			inputs(data);
			inputs<Other...>(other...);
			return *this;
		}

		template <class TParameter>
		Link& outputs(TParameter &data)
		{
			if (data.version[interval] < version) {
				data.version[interval] = version;
				++update;
			}
			data.isset = isset;
			return *this;
		}

		template <class TParameter, class ...Other>
		Link& outputs(TParameter &data, Other&...other)
		{
			outputs(data);
			outputs<Other...>(other...);
			return *this;
		}
	};

	const int interval, isconst, update;

	Operator(int interval, bool isconst, bool update): interval(interval), isconst(isconst), update(update) {}
	virtual ~Operator() {}

	virtual const char* name() =0;
	virtual void operator++() =0;
};

struct OperatorBuilder {
	virtual ~OperatorBuilder() {}

	virtual const char* name() =0;
	virtual void now() =0;

	virtual bool build(HeatTransferKernelOpt &kernel) { return true; }
};

struct ElementOperatorBuilder: public OperatorBuilder {
	virtual ~ElementOperatorBuilder() {}

	void now();

	virtual void apply(int interval) =0;
};

struct BoundaryOperatorBuilder: public OperatorBuilder {
	virtual ~BoundaryOperatorBuilder() {}

	void now();

	virtual void apply(int region, int interval) =0;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */