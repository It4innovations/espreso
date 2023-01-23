
#include "w.papi.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"

#include <cstdio>
#include <vector>

using namespace espreso;

int PAPI::status = 0;
int PAPI::values = 0;

#ifdef HAVE_PAPI
#include "papi.h"

PAPI::PAPI(): set{PAPI_NULL}
{

}

PAPI::~PAPI()
{
	std::vector<long long> dummy(values);
	PAPI_stop(set, dummy.data());
	PAPI_cleanup_eventset(set);
	PAPI_destroy_eventset(&set);
}

void PAPI::init()
{
	status += PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT;
	status += PAPI_create_eventset(&set);
	if (info::ecf->output.papi_events.size() && status == 0) {
		std::vector<std::string> events = Parser::split(info::ecf->output.papi_events, " ");
		for (size_t i = 0; status == 0 && i < events.size(); ++i) {
			int native;
			status += PAPI_event_name_to_code(const_cast<char*>(events[i].c_str()), &native);
			status += PAPI_add_event(set, native);
		}
		values = events.size();
	}
	if (info::ecf->output.papi_codes.size() && status == 0) {
		std::vector<std::string> codes = Parser::split(info::ecf->output.papi_codes, " ");
		for (size_t i = 0; status == 0 && i < codes.size(); ++i) {
			int code = strtol(codes[i].c_str(), NULL, 16);
			status += PAPI_add_event(set, code);
		}
		values = codes.size();
	}
	status += PAPI_start(set);
}

void PAPI::read(long long *values)
{
	if (set != PAPI_NULL) {
		PAPI_read(set, values);
	}
}

#else

PAPI::PAPI(): set{0} {}
PAPI::~PAPI() {}
void PAPI::init() {}
void PAPI::read(long long *values) {}

#endif
