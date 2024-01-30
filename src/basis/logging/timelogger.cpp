
#include "timelogger.h"
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

double TimeLogger::initClockTime = TimeLogger::time();
time_t TimeLogger::initTime = std::time(NULL);

double TimeLogger::time()
{
	return omp_get_wtime();
}

double TimeLogger::duration()
{
	return omp_get_wtime() - initClockTime;
}

static void mergeEvents(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		TimeLogger::Event::Data &inmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->min;
		TimeLogger::Event::Data &inmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->max;
		TimeLogger::Event::Data &inavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->avg;
		TimeLogger::Event::Data &indmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->dmin;
		TimeLogger::Event::Data &indmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->dmax;
		TimeLogger::Event::Data &indavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->davg;
		TimeLogger::Event::Data &insmin = (static_cast<TimeLogger::EventStatistics*>(in) + i)->smin;
		TimeLogger::Event::Data &insmax = (static_cast<TimeLogger::EventStatistics*>(in) + i)->smax;
		TimeLogger::Event::Data &insavg = (static_cast<TimeLogger::EventStatistics*>(in) + i)->savg;
		TimeLogger::Event::Data &outmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->min;
		TimeLogger::Event::Data &outmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->max;
		TimeLogger::Event::Data &outavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->avg;
		TimeLogger::Event::Data &outdmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->dmin;
		TimeLogger::Event::Data &outdmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->dmax;
		TimeLogger::Event::Data &outdavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->davg;
		TimeLogger::Event::Data &outsmin = (static_cast<TimeLogger::EventStatistics*>(out) + i)->smin;
		TimeLogger::Event::Data &outsmax = (static_cast<TimeLogger::EventStatistics*>(out) + i)->smax;
		TimeLogger::Event::Data &outsavg = (static_cast<TimeLogger::EventStatistics*>(out) + i)->savg;

		switch ((static_cast<TimeLogger::EventStatistics*>(out) + i)->type) {
		case TimeLogger::Event::START:
		case TimeLogger::Event::CHECKPOINT:
		case TimeLogger::Event::ACCUMULATED:
		case TimeLogger::Event::END:
			outmin.time = std::min(outmin.time, inmin.time);
			outmax.time = std::max(outmax.time, inmax.time);
			outavg.time += inavg.time;
			outdmin.time = std::min(outdmin.time, indmin.time);
			outdmax.time = std::max(outdmax.time, indmax.time);
			outdavg.time += indavg.time;
			outsmin.time = std::min(outsmin.time, insmin.time);
			outsmax.time = std::max(outsmax.time, insmax.time);
			outsavg.time += insavg.time;
			break;
		case TimeLogger::Event::INT:
			outmin.ivalue = std::min(outmin.ivalue, inmin.ivalue);
			outmax.ivalue = std::max(outmax.ivalue, inmax.ivalue);
			outavg.ivalue += inavg.ivalue;
			break;
		case TimeLogger::Event::LONG:
			outmin.lvalue = std::min(outmin.lvalue, inmin.lvalue);
			outmax.lvalue = std::max(outmax.lvalue, inmax.lvalue);
			outavg.lvalue += inavg.lvalue;
			break;
		case TimeLogger::Event::SIZE:
			outmin.svalue = std::min(outmin.svalue, inmin.svalue);
			outmax.svalue = std::max(outmax.svalue, inmax.svalue);
			outavg.svalue += inavg.svalue;
			break;
		case TimeLogger::Event::DOUBLE:
			outmin.dvalue = std::min(outmin.dvalue, inmin.dvalue);
			outmax.dvalue = std::max(outmax.dvalue, inmax.dvalue);
			outavg.dvalue += inavg.dvalue;
			break;
		default:
			break;
		}
	}
}

static void printdata(
		const char* format, const char* name, const char* suffix,
		double avg, double min, double max, double sectiontime)
{
	std::string fullname = std::string(name) + std::string(suffix);
	avg /= info::mpi::size;
	sectiontime /= info::mpi::size;
	std::string savg = std::to_string(avg);
	std::string smin = std::to_string(min);
	std::string smax = std::to_string(max);
	char sratio[6], simb[6];
	snprintf(sratio, 6, "%5.2f", 100 * avg / sectiontime);
	snprintf(simb, 6, "%5.2f", max / min);
	savg[8] = smin[8] = smax[8] = '\0';
	sratio[5] = simb[5] = '\0';
	eslog::info(format, fullname.c_str(), savg.c_str(), smin.c_str(), smax.c_str(), sratio, simb);
}

void TimeLogger::finish()
{
	if (verbosity == 0) {
		return;
	}
	profiler::syncstart("time_statistics");

	std::vector<EventStatistics> events(_events.begin(), _events.end());
	std::vector<double> prev(10), begin(10);
	std::vector<std::vector<const char*> > accumulated;
	std::vector<size_t> block;

	size_t namewidth = 43;

	for (size_t i = 0; i < events.size(); i++) {
		switch (events[i].type) {
		case Event::START:
			block.push_back(i);
			prev.push_back(events[i].data.time);
			begin.push_back(events[i].data.time);
			events[i].data.time -= initClockTime;
			events[i].duration.time = 0;
			accumulated.push_back({});
			break;
		case Event::ACCUMULATED:
			if (std::find(accumulated.back().begin(), accumulated.back().end(), events[i].name) == accumulated.back().end()) {
				accumulated.back().push_back(events[i].name);
			}
			/* no break */
		case Event::CHECKPOINT:
			events[i].data.time -= prev.back();
			events[i].duration.time -= begin.back();
			events[i].sum.time = events[i].data.time;
			prev.back() += events[i].data.time;
			namewidth = std::max(namewidth, strlen(events[i].name));
			break;
		case Event::END:
			events[i].data.time -= prev.back();
			events[i].duration.time -= begin.back();
			prev.pop_back();
			begin.pop_back();
			namewidth = std::max(namewidth, strlen(events[i].name));
			for (auto acc = accumulated.back().begin(); acc != accumulated.back().end(); ++acc) {
				double sum = 0;
				for (size_t j = block.back(); j < i; ++j) {
					if (*acc == events[j].name) {
						sum += events[j].data.time;
						// it is enough to have final sum in the last accumulated event
						events[j].savg.time = events[j].smin.time = events[j].smax.time = events[j].sum.time = sum;
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
		events[i].min = events[i].data;
		events[i].max = events[i].data;
		events[i].avg = events[i].data;
		events[i].dmin = events[i].duration;
		events[i].dmax = events[i].duration;
		events[i].davg = events[i].duration;
	}

	std::vector<EventStatistics> statistics(events);

	{ // synchronize across processes
		size_t eventsize = events.size(), minsize, maxsize;
		Communication::allReduce(&eventsize, &minsize, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
		Communication::allReduce(&eventsize, &maxsize, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
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
		profiler::syncend("time_statistics");
		return;
	}

	double duration = TimeLogger::duration() * info::mpi::size;
	const char* usershead  = "  %-43s %s  <%s - %s> [%s] [%s]\n";
	const char* usersdata  = "   %-42s %s  <%s - %s> [%s] [%s]\n";
	const char* parserformat = "time: region='%-42s', avg='%s', min='%s', max='%s', ratio='%s', imb='%s'\n";
	const char *headformat = NULL, *dataformat = NULL;

	switch (info::ecf->output.logger) {
	case OutputConfiguration::LOGGER::USER: headformat = usershead; dataformat = usersdata; break;
	case OutputConfiguration::LOGGER::PARSER: headformat = parserformat; dataformat = parserformat; break;
	}

	auto print = [&] (size_t start, size_t end, int printeddepth) {
		eslog::info(" ============================================ avg. [s]  < min [s] -  max [s]> [  %  ] [ imb ]   \n");
		int depth = printeddepth - 1, loadstep = 1;
		int allowparam = 0, denyparam = 1;
		for (size_t i = start; i <= end; i++) {
			switch (statistics[i].type) {
			case Event::START:
				denyparam = 0;
				if (allowparam) {
					printdata( headformat, events[start].name, " TOTAL DURATION",
							statistics[allowparam].davg.time, statistics[allowparam].dmin.time, statistics[allowparam].dmax.time, duration);
					allowparam = 0;
					--depth;
				}
				++depth;
				if (depth == printeddepth) {
					printdata(headformat, events[i].name, " STARTED AT",
							statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, duration);
				}
				break;
			case Event::CHECKPOINT:
				denyparam = 0;
				if (allowparam) {
					printdata( headformat, events[start].name, " TOTAL DURATION",
							statistics[allowparam].davg.time, statistics[allowparam].dmin.time, statistics[allowparam].dmax.time, duration);
					allowparam = 0;
					--depth;
				}
				if (depth == printeddepth) {
					printdata(dataformat, events[i].name, "",
							statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].davg.time);
				}
				break;
			case Event::ACCUMULATED:
				break;
			case Event::END:
				denyparam = 0;
				if (depth == printeddepth) {
					printdata( dataformat, events[i].name, "",
							statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].davg.time);
					allowparam = i;
				} else {
					denyparam = i;
				}
				if (allowparam == 0) {
					--depth;
				}
				break;
			case Event::INT:
				if (depth > 1 && ((depth == printeddepth && !denyparam) || allowparam)) {
					eslog::info("    [param=%s] %*d  <%8d - %8d>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.ivalue / info::mpi::size,
							statistics[i].min.ivalue, statistics[i].max.ivalue,
							(double)statistics[i].max.ivalue / statistics[i].min.ivalue);
				}
				break;
			case Event::LONG:
				if (depth > 1 && ((depth == printeddepth && !denyparam) || allowparam)) {
					eslog::info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.lvalue / info::mpi::size,
							statistics[i].min.lvalue, statistics[i].max.lvalue,
							(double)statistics[i].max.lvalue / statistics[i].min.lvalue);
				}
				break;
			case Event::SIZE:
				if (depth > 1 && ((depth == printeddepth && !denyparam) || allowparam)) {
					eslog::info("    [param=%s] %*ld  <%8ld - %8ld>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							statistics[i].avg.svalue / info::mpi::size,
							statistics[i].min.svalue, statistics[i].max.svalue,
							(double)statistics[i].max.svalue / statistics[i].min.svalue);
				}
				break;
			case Event::DOUBLE:
				if (depth > 1 && ((depth == printeddepth && !denyparam) || allowparam)) {
					char avg[9], min[9], max[9];
					snprintf(avg, 9, "%f", statistics[i].avg.dvalue / info::mpi::size);
					snprintf(min, 9, "%f", statistics[i].min.dvalue);
					snprintf(max, 9, "%f", statistics[i].max.dvalue);
					eslog::info("    [param=%s] %*s  <%s - %s>         [%5.2f]\n",
							events[i].name, namewidth - strlen(events[i].name) - 2,
							avg, min, max,
							statistics[i].max.dvalue / statistics[i].min.dvalue);
				}
				break;
			case Event::LOADSTEP:
				if (depth == 0) {
					eslog::info("                                                                              [ LOAD STEP%2d ]\n", loadstep++);
				}
				denyparam = 0;
				if (allowparam) {
					printdata( headformat, events[start].name, " TOTAL DURATION",
							statistics[allowparam].davg.time, statistics[allowparam].dmin.time, statistics[allowparam].dmax.time, duration);
					allowparam = 0;
					--depth;
				}
				break;
			default:
				break;
			}
		}
		if (allowparam) {
			printdata( headformat, events[start].name, " TOTAL DURATION",
					statistics[allowparam].davg.time, statistics[allowparam].dmin.time, statistics[allowparam].dmax.time, duration);
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
							double min = duration, max = 0, savg = 0, smin = 0, smax = 0;
							for (size_t j = start; j <= end; j++) {
								if (events[i].name == events[j].name) {
									min = std::min(min, statistics[j].min.time);
									max = std::max(max, statistics[j].max.time);
									savg = statistics[j].savg.time;
									smin = statistics[j].smin.time;
									smax = statistics[j].smax.time;
									++counter;
								}
							}
							eslog::info("   %s [%dx]\n", events[i].name, counter);
							printdata(dataformat, events[i].name, " [1]", statistics[i].avg.time, statistics[i].min.time, statistics[i].max.time, statistics[end].davg.time);
							printdata(dataformat, events[i].name, " [~]", savg / counter, min, max, statistics[end].davg.time);
							printdata(dataformat, events[i].name, " [+]", savg, smin, smax, statistics[end].davg.time);
						}
					}
				}
				break;
			case Event::END:
				--depth;
				break;
			case Event::INT:
			case Event::LONG:
			case Event::SIZE:
			case Event::DOUBLE:
			case Event::LOADSTEP:
			default:
				break;
			}
		}
		eslog::info(" ============================================================================================= \n");
	};

	// WARNING: MPI sometimes copy name pointer from other process, hence use _events names
	std::vector<size_t> begins;
	std::vector<std::vector<const char*> > uniques, duplications;
	int loadstep = 0;
	size_t lastend = 0;

	auto isparam = [&] (size_t i) {
		return
				statistics[i].type == Event::INT ||
				statistics[i].type == Event::LONG ||
				statistics[i].type == Event::SIZE ||
				statistics[i].type == Event::DOUBLE;
	};

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
			while (i + 1 < statistics.size() && isparam(i + 1)) {
				++i;
			}
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

	eslog::info(" == OVERALL TIME ============================================================================= \n");
	print(0, lastend, 0);
	Communication::barrier();
	profiler::syncend("time_statistics");
}



