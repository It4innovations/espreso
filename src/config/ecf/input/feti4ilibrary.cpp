
#include "feti4ilibrary.h"
#include "config/configuration.hpp"

espreso::FETI4ILibraryConfiguration::FETI4ILibraryConfiguration()
{
	domains = 4;
	REGISTER(domains, ECFMetaData()
			.setdescription({ "Number of domains in each cluster." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(solver, ECFMetaData()
			.setdescription({ "Settings of FETI4I linear solver." }));

	for (size_t i = 0; i < solver.ecfdescription->parameters.size(); ) {
		if (
				solver.ecfdescription->parameters[i]->data() == &solver.max_iterations ||
				solver.ecfdescription->parameters[i]->data() == &solver.method ||
				solver.ecfdescription->parameters[i]->data() == &solver.preconditioner ||
				solver.ecfdescription->parameters[i]->data() == &solver.iterative_solver ||
				solver.ecfdescription->parameters[i]->data() == &solver.precision ||
				solver.ecfdescription->parameters[i]->data() == &solver.sc_size) {

			// keep only the above
			i++;
		} else {
			solver.ecfdescription->dropParameter(solver.ecfdescription->parameters[i]);
		}
	}
}



