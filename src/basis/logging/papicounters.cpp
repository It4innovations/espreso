
#include "papicounters.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"
#include "wrappers/mpi/communication.h"

#include "omp.h"
#include <cstdio>
#include <cstring>
#include <ctime>

using namespace espreso;

static void mergeEvents(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		PAPICounters::EventStatistics &ein = *(static_cast<PAPICounters::EventStatistics*>(in) + i);
		PAPICounters::EventStatistics &eout = *(static_cast<PAPICounters::EventStatistics*>(out) + i);

		switch ((static_cast<PAPICounters::EventStatistics*>(out) + i)->type) {
		case PAPICounters::Event::START:
		case PAPICounters::Event::CHECKPOINT:
		case PAPICounters::Event::ACCUMULATED:
		case PAPICounters::Event::END:
			eout.min = std::min(eout.min, ein.min);
			eout.max = std::max(eout.max, ein.max);
			eout.avg += ein.avg;
			eout.dmin = std::min(eout.dmin, ein.dmin);
			eout.dmax = std::max(eout.dmax, ein.dmax);
			eout.davg += ein.davg;
			eout.smin = std::min(eout.smin, ein.smin);
			eout.smax = std::max(eout.smax, ein.smax);
			eout.savg += ein.savg;
			break;
		default:
			break;
		}
	}
}

static void printdata(
		const char* format, const char* name, const char* suffix,
		long avg, long min, long max, long sectionTotal)
{
	std::string fullname = std::string(name) + std::string(suffix);
	char sratio[6] = "   ", simb[6] = "   ";
	snprintf(sratio, 5, "%3d", (int)(100.0 * avg / sectionTotal));
	if (min != 0) {
		snprintf(simb, 5, "%3d", (int)(max / min));
	}
	sratio[5] = simb[5] = '\0';
	eslog::info(format, fullname.c_str(), (double)avg / info::mpi::size, min, max, sratio, simb);
}

void PAPICounters::finish()
{
	if (verbosity == 0 || (info::ecf->output.papi_code == 0 && info::ecf->output.papi_event.size() == 0)) {
		return;
	}
	if (!papi.isValid()) {
		eslog::info("\n ============================================================================================= \n");
		if (info::ecf->output.papi_event.size()) {
			eslog::info(" == INVALID PAPI EVENT %68s == \n", info::ecf->output.papi_event.c_str());
		}
		if (info::ecf->output.papi_code) {
			eslog::info(" == INVALIED PAPI CODE                                                           0x%8x == \n", info::ecf->output.papi_code);
		}
		eslog::info(" ============================================================================================= \n");
		return;
	}
	long duration = papi.read() - init;
	std::vector<EventStatistics> events(_events.begin(), _events.end());
	std::vector<long> prev(10), begin(10);
	std::vector<std::vector<const char*> > accumulated;
	std::vector<size_t> block;

	size_t namewidth = 43;

	for (size_t i = 0; i < events.size(); i++) {
		switch (events[i].type) {
		case Event::START:
			block.push_back(i);
			prev.push_back(events[i].value);
			begin.push_back(events[i].value);
			events[i].value -= init;
			events[i].duration = 0;
			accumulated.push_back({});
			break;
		case Event::ACCUMULATED:
			if (std::find(accumulated.back().begin(), accumulated.back().end(), events[i].name) == accumulated.back().end()) {
				accumulated.back().push_back(events[i].name);
			}
			/* no break */
		case Event::CHECKPOINT:
			events[i].value -= prev.back();
			events[i].duration -= begin.back();
			events[i].sum = events[i].value;
			prev.back() += events[i].value;
			namewidth = std::max(namewidth, strlen(events[i].name));
			break;
		case Event::END:
			events[i].value -= prev.back();
			events[i].duration -= begin.back();
			prev.pop_back();
			begin.pop_back();
			namewidth = std::max(namewidth, strlen(events[i].name));
			for (auto acc = accumulated.back().begin(); acc != accumulated.back().end(); ++acc) {
				double sum = 0;
				for (size_t j = block.back(); j < i; ++j) {
					if (*acc == events[j].name) {
						sum += events[j].value;
						// it is enough to have final sum in the last accumulated event
						events[j].sum = sum;
					}
				}
			}
			accumulated.pop_back();
			block.pop_back();
			break;
		default:
			namewidth = std::max(namewidth, strlen(events[i].name) + 2);
			break;
		}
		events[i].min = events[i].value;
		events[i].max = events[i].value;
		events[i].avg = events[i].value;
		events[i].dmin = events[i].duration;
		events[i].dmax = events[i].duration;
		events[i].davg = events[i].duration;
	}

	std::vector<EventStatistics> statistics(events);

	{ // synchronize across processes
		size_t eventsize = events.size(), minsize, maxsize;
		Communication::allReduce(&eventsize, &minsize, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
		Communication::allReduce(&eventsize, &maxsize, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
		Communication::allReduce(&duration, nullptr, 1, MPITools::getType(duration).mpitype, MPI_SUM);
		if (minsize == eventsize && maxsize == eventsize) {
			MPI_Op reduce;
			MPI_Datatype mpievent;
			MPI_Op_create(mergeEvents, 1, &reduce);
			MPI_Type_contiguous(sizeof(EventStatistics), MPI_BYTE, &mpievent);
			MPI_Type_commit(&mpievent);
			Communication::reduce(events.data(), statistics.data(), events.size(), mpievent, reduce, 0);
			MPI_Op_free(&reduce);
			MPI_Type_free(&mpievent);
		} else {
			eslog::warning("Various number of time events (only root data are printed)\n");
		}
	}

	if (grank) {
		Communication::barrier();
		return;
	}

	const char* usershead  = "  %-41s %9.3g %13lu %13lu [%s] [%s]\n";
	const char* usersdata  = "   %-40s %9.3g %13lu %13lu [%s] [%s]\n";
	const char* parserformat = "papi: region='%-41s', avg='%g', min='%lu', max='%lu', ratio='%s', imb='%s'\n";
	const char *headformat = NULL, *dataformat = NULL;

	switch (info::ecf->output.logger) {
	case OutputConfiguration::LOGGER::USER: headformat = usershead; dataformat = usersdata; break;
	case OutputConfiguration::LOGGER::PARSER: headformat = parserformat; dataformat = parserformat; break;
	}

	auto print = [&] (size_t start, size_t end, int printeddepth) {
		eslog::info(" ==========================================       avg           min           max [ % ] [imb]   \n");
		int depth = printeddepth - 1, loadstep = 1;
		int allowparam = 0;
		for (size_t i = start; i <= end; i++) {
			switch (statistics[i].type) {
			case Event::START:
				if (allowparam) {
					printdata(headformat, events[start].name, " TOTAL", statistics[allowparam].davg, statistics[allowparam].dmin, statistics[allowparam].dmax, duration);
					allowparam = 0;
					--depth;
				}
				++depth;
				if (depth == printeddepth) {
					printdata(headformat, events[i].name, " STARTED AT", statistics[i].avg, statistics[i].min, statistics[i].max, duration);
				}
				break;
			case Event::CHECKPOINT:
				if (allowparam) {
					printdata(headformat, events[start].name, " TOTAL", statistics[allowparam].davg, statistics[allowparam].dmin, statistics[allowparam].dmax, duration);
					allowparam = 0;
					--depth;
				}
				if (depth == printeddepth) {
					printdata(dataformat, events[i].name, "", statistics[i].avg, statistics[i].min, statistics[i].max, statistics[end].davg);
				}
				break;
			case Event::ACCUMULATED:
				break;
			case Event::END:
				if (depth == printeddepth) {
					printdata(dataformat, events[i].name, "", statistics[i].avg, statistics[i].min, statistics[i].max, statistics[end].davg);
					allowparam = i;
				}
				if (allowparam == 0) {
					--depth;
				}
				break;
			case Event::LOADSTEP:
				if (depth == 0) {
					eslog::info("                                                                              [ LOAD STEP%2d ]\n", loadstep++);
				}
				if (allowparam) {
					printdata(headformat, events[start].name, " TOTAL", statistics[allowparam].davg, statistics[allowparam].dmin, statistics[allowparam].dmax, duration);
					allowparam = 0;
					--depth;
				}
				break;
			default:
				break;
			}
		}
		if (allowparam) {
			printdata(headformat, events[start].name, " TOTAL", statistics[allowparam].davg, statistics[allowparam].dmin, statistics[allowparam].dmax, duration);
		}
		eslog::info(" ============================================================================================= \n");
	};

	auto printduplications = [&] (size_t start, size_t end, int printeddepth, const std::vector<const char*> &duplications) {
		int depth = printeddepth - 1;
		std::vector<const char*> printed;
		for (size_t i = start; i <= end; i++) {
			switch (statistics[i].type) {
			case Event::START:
				++depth;
				break;
			case Event::CHECKPOINT:
				break;
			case Event::ACCUMULATED:
				if (depth == printeddepth) {
					if (std::find(duplications.begin(), duplications.end(), events[i].name) != duplications.end()) {
						if (std::find(printed.begin(), printed.end(), events[i].name) == printed.end()) {
							printed.push_back(events[i].name);
							int counter = 0;
							long min = duration, max = 0, savg, smin, smax;
							for (size_t j = start; j <= end; j++) {
								if (events[i].name == events[j].name) {
									min = std::min(min, statistics[j].min);
									max = std::max(max, statistics[j].max);
									savg = statistics[j].savg;
									smin = statistics[j].smin;
									smax = statistics[j].smax;
									++counter;
								}
							}
							eslog::info("   %s [%dx]\n", events[i].name, counter);
							printdata(dataformat, events[i].name, " [1]", statistics[i].avg, statistics[i].min, statistics[i].max, statistics[end].avg);
							printdata(dataformat, events[i].name, " [~]", savg / counter, min, max, statistics[end].avg);
							printdata(dataformat, events[i].name, " [+]", savg, smin, smax, statistics[end].avg);
						}
					}
				}
				break;
			case Event::END:
				--depth;
				break;
			default:
				break;
			}
		}
		eslog::info(" ============================================================================================= \n");
	};

	// WARNING: MPI sometimes copy name pointer from other process, hence use _events names
	std::vector<size_t> begins;
	std::vector<std::vector<const char*> > uniques, duplications;
	std::vector<std::string> solvers(step::step.loadstep + 1);
	int loadstep = 0;
	size_t lastend = 0;

	for (size_t i = 0; i < statistics.size(); i++) {
		switch (statistics[i].type) {
		case Event::START:
			begins.push_back(i);
			uniques.push_back({});
			duplications.push_back({});
			break;
		case Event::ACCUMULATED:
		case Event::CHECKPOINT:
			if (std::find(uniques.back().begin(), uniques.back().end(), events[i].name) == uniques.back().end()) {
				uniques.back().push_back(events[i].name);
			} else {
				if (std::find(duplications.back().begin(), duplications.back().end(), events[i].name) == duplications.back().end()) {
					duplications.back().push_back(events[i].name);
				}
			}
			break;
		case Event::END:
			lastend = i;
			if (begins.size() > 1) {
				print(begins.back(), i, begins.size() + 1);
				printduplications(begins.back(), i, begins.size() + 1, duplications.back()); // do not print loadsteps
			}
			begins.pop_back();
			uniques.pop_back();
			duplications.pop_back();
			break;
		case Event::LOADSTEP:
			++loadstep;
			if (verbosity > 1) {
				eslog::info(" == LOADSTEP: %2d ============================================================================= \n", loadstep);
			}
			break;
		default:
			break;
		}
	}

	eslog::info("\n ============================================================================================= \n");
	if (info::ecf->output.papi_event.size() && PAPI::isValid()) {
		eslog::info(" == PAPI EVENT %76s == \n", info::ecf->output.papi_event.c_str());
	}
	if (info::ecf->output.papi_code && PAPI::isValid()) {
		eslog::info(" == PAPI CODE                                                                    0x%8x == \n", info::ecf->output.papi_code);
	}
	eslog::info(" ============================================================================================= \n");
	print(0, lastend, 0);
	Communication::barrier();
	profiler::syncend("time_statistics");
}
