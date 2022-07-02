
#include "w.backward.h"
#include "esinfo/eslog.hpp"

#include <string>
#include <sstream>

#ifdef HAVE_BACKWARD
#include "backward.hpp"
#endif

void espreso::BackWard::call()
{
#ifdef HAVE_BACKWARD
	backward::StackTrace st;
	st.load_here(32);
	std::stringstream ss;
	backward::Printer p;
	p.snippet = true;
	p.color_mode = backward::ColorMode::type::always;
	p.address = true;
	p.object = false;
	p.print(st, ss);

	while (ss.good()) {
		std::string line;
		getline(ss, line);
		eslog::info("%s\n", line.c_str());
	}
#endif
}
