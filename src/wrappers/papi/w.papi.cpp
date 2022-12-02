
#include "w.papi.h"
#include "esinfo/ecfinfo.h"

#include <cstdio>

using namespace espreso;

int PAPI::status = 0;

#ifdef HAVE_PAPI
#include "papi.h"

PAPI::PAPI(): set{PAPI_NULL}
{

}

PAPI::~PAPI()
{
	long long value;
	PAPI_stop(set, &value);
	PAPI_cleanup_eventset(set);
	PAPI_destroy_eventset(&set);
}

void PAPI::init()
{
	status += PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT;
	status += PAPI_create_eventset(&set);
	if (info::ecf->output.papi_event.size() && PAPI::isValid()) {
		int native;
		status += PAPI_event_name_to_code(const_cast<char*>(info::ecf->output.papi_event.c_str()), &native);
		status += PAPI_add_event(set, native);
	}
	if (info::ecf->output.papi_code && PAPI::isValid()) {
		status += PAPI_add_event(set, info::ecf->output.papi_code);
	}
	status += PAPI_start(set);
}

long PAPI::read()
{
	long long value = 0;
	if (set != PAPI_NULL) {
		PAPI_read(set, &value);
	}
	return value;
}

#else

PAPI::PAPI(): set{0} {}
PAPI::~PAPI() {}
void PAPI::init(int event) {}
long long PAPI::read() { return 0; }

#endif
