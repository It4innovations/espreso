
#include <config/ecf/ecf.h>
#include "output.h"
#include "config/configuration.hpp"

using namespace espreso;

MonitorConfiguration::MonitorConfiguration(const ECF *root)
: _root(root)
{
	REGISTER(region, ECFMetaData()
			.setdescription({ "Region" })
			.setdatatype({ ECFDataType::REGION }));

	statistics = STATISTICS::AVG;
	REGISTER(statistics, ECFMetaData()
			.setdescription({ "Statistics" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("MIN").setdescription("Minimum."))
			.addoption(ECFOption().setname("MAX").setdescription("Maximum."))
			.addoption(ECFOption().setname("AVG").setdescription("Average."))
			.addoption(ECFOption().setname("NORM").setdescription("Norm (for testing purposes).").allowonly([] () { return false; }))
			.addoption(ECFOption().setname("ABSMIN").setdescription("Minimum of absolute value."))
			.addoption(ECFOption().setname("ABSMAX").setdescription("Maximum of absolute value.")));

	REGISTER(property, ECFMetaData()
			.setdescription({ "Result" })
			.setdatatype({ ECFDataType::OPTION }));

	HeatTransferOutputSettings::addMonitorableProperties(ecfdescription->getParameter(&property)->metadata, _root);
	StructuralMechanicsOutputSettings::addMonitorableProperties(ecfdescription->getParameter(&property)->metadata, _root);
}

ResultsSelectionConfiguration::ResultsSelectionConfiguration(const PhysicsConfiguration::TYPE &physics)
: _physics(physics)
{
	basic();

	REGISTER(thickness, ECFMetaData()
			.setdescription({ "Element thickness." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D || _physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D; }));
}

HarmonicOuputConfiguration::HarmonicOuputConfiguration()
{
	samples = 20;
	REGISTER(samples, ECFMetaData()
			.setdescription({ "Number of samples for each frequency." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	results_store_frequency = STORE_FREQUENCY::NEVER;
	results_nth_stepping = 10;
	REGISTER(results_store_frequency, ECFMetaData()
			.setdescription({ "Results store frequency" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEVER").setdescription("Storing is turned off."))
			.addoption(ECFOption().setname("EVERY_FREQUENCY").setdescription("Results are stored for each frequency."))
			.addoption(ECFOption().setname("EVERY_NTH_FREQUENCY").setdescription("Results are stored for each nth frequency."))
			.addoption(ECFOption().setname("SPECIFIC_FREQUENCIES").setdescription("Only requested frequencies are stored.")));

	REGISTER(results_nth_stepping, ECFMetaData()
			.setdescription({ "Results store n-th stepping." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(requested_frequencies, ECFMetaData()
			.setdescription({ "Requested frequencies." })
			.setdatatype({ ECFDataType::STRING }));
}

OutputConfiguration::OutputConfiguration(const ECF *root)
: results_selection(root->physics), _root(root)
{
	log_dir = "debug";
	REGISTER(log_dir, ECFMetaData()
			.setdescription({ "A name of logging directory" })
			.setdatatype({ ECFDataType::STRING }));

	verbose_level = 1;
	REGISTER(verbose_level, ECFMetaData()
			.setdescription({ "Verbose level [0-3]." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	measure_level = 1;
	REGISTER(measure_level, ECFMetaData()
			.setdescription({ "Measure level [0-3]." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(papi_codes, ECFMetaData()
			.setdescription({ "List of reported PAPI event codes." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(papi_events, ECFMetaData()
			.setdescription({ "List of reported PAPI events." })
			.setdatatype({ ECFDataType::STRING }));

	logger = LOGGER::USER;
	REGISTER(logger, ECFMetaData()
			.setdescription({ "Outpu logger settings" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("USER").setdescription("Output for users."))
			.addoption(ECFOption().setname("PARSER").setdescription("A parser frienly output.")));

	print_matrices = 0;
	REGISTER(print_matrices, ECFMetaData()
			.setdescription({ "Print assembler matrices for debugging." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	store_decomposition = 1;
	REGISTER(store_decomposition, ECFMetaData()
			.setdescription({ "Store decomposition." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	path = "results";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path" })
			.setdatatype({ ECFDataType::STRING }));

	format = FORMAT::ENSIGHT;
	REGISTER(format, ECFMetaData()
			.setdescription({ "Format" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VTK_LEGACY").setdescription("VTK legacy format."))
//			.addoption(ECFOption().setname("VTK_XML_ASCII").setdescription("VTK XML ASCII format."))
//			.addoption(ECFOption().setname("VTK_XML_BINARY").setdescription("VTK XML binary format."))
			.addoption(ECFOption().setname("ENSIGHT").setdescription("EnSight format."))
			.addoption(ECFOption().setname("XDMF").setdescription("XDMF format."))
			.addoption(ECFOption().setname("STL_SURFACE").setdescription("Surface of bodies in STL format."))
			.addoption(ECFOption().setname("NETGEN").setdescription("Netgen neutral format (only for tetrahedral meshes).")));

	mode = MODE::PTHREAD;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "ASYNC library mode" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SYNC").setdescription("Output is synchronized."))
			.addoption(ECFOption().setname("PTHREAD").setdescription("Output is done by asynchronous thread.")));

	writer = WRITER::MPI;
	REGISTER(writer, ECFMetaData()
			.setdescription({ "Set the type of writter." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("MPI").setdescription("Output is write by MPI_file_write_at."))
			.addoption(ECFOption().setname("MPI_COLLECTIVE").setdescription("Output is write by MPI_file_write_at_all."))
//			.addoption(ECFOption().setname("POSIX").setdescription("Output is write by POSIX API.")) // POSIX does not work
			);

	stripe_size = 1024 * 1024; // 1M
	REGISTER(stripe_size, ECFMetaData()
				.setdescription({ "The size of data to be stored by each writter." })
				.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	stripe_count = 0; // all processes store its own parts
	REGISTER(stripe_count, ECFMetaData()
				.setdescription({ "The number of writters." })
				.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

//	output_node_group_size = 7;
//	REGISTER(output_node_group_size, ECFMetaData()
//			.setdescription({ "Number of MPI processes that send output data to the same storing node." })
//			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSeparator();

	results_store_frequency = monitors_store_frequency = STORE_FREQUENCY::EVERY_SUBSTEP;
	results_nth_stepping = monitors_nth_stepping = 10;
	REGISTER(results_store_frequency, ECFMetaData()
			.setdescription({ "Results store frequency" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEVER").setdescription("Storing is turned off."))
			.addoption(ECFOption().setname("EVERY_SUBSTEP").setdescription("Results are stored after each time step."))
			.addoption(ECFOption().setname("EVERY_NTH_SUBSTEP").setdescription("Results are stored after each nth time step."))
			.addoption(ECFOption().setname("LAST_SUBSTEP").setdescription("Only last results are stored."))
			.addoption(ECFOption().setname("DEBUG").setdescription("Storing also iteration results.")));
	REGISTER(monitors_store_frequency, ECFMetaData()
			.setdescription({ "Monitoring store frequency" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEVER").setdescription("Storing is turned off."))
			.addoption(ECFOption().setname("EVERY_SUBSTEP").setdescription("Monitors are stored after each time step."))
			.addoption(ECFOption().setname("EVERY_NTH_SUBSTEP").setdescription("Monitors are stored after each nth time step."))
			.addoption(ECFOption().setname("LAST_SUBSTEP").setdescription("Only last monitors are stored.")));

	REGISTER(results_nth_stepping, ECFMetaData()
			.setdescription({ "Write results" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.allowonly([&] () { return results_store_frequency == STORE_FREQUENCY::EVERY_NTH_SUBSTEP; }));
	REGISTER(monitors_nth_stepping, ECFMetaData()
			.setdescription({ "Monitors store stepping" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER})
			.allowonly([&] () { return monitors_store_frequency == STORE_FREQUENCY::EVERY_NTH_SUBSTEP; }));
	REGISTER(frequency_to_time, ECFMetaData()
			.setdescription({ "Frequency to time settings" })
			.allowonly([&] () { return _root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }));

	ecfdescription->addSpace();


	store_results = STORE_RESULTS::BASIC;
	REGISTER(store_results, ECFMetaData()
			.setdescription({ "Stored properties" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BASIC").setdescription("Basic properties."))
			.addoption(ECFOption().setname("ALL").setdescription("All properties."))
			.addoption(ECFOption().setname("USER").setdescription("User defined properties.")))
	->addListener(ECFParameter::Event::VALUE_SET, [&] (const std::string &value) {
		switch (store_results) {
		case STORE_RESULTS::BASIC:
			results_selection.basic();
			break;
		case STORE_RESULTS::ALL:
			results_selection.all();
			break;
		case STORE_RESULTS::USER:
			break;
		}
	});

	REGISTER(results_selection, ECFMetaData()
			.setdescription({ "Properties selection" })
			.allowonly([&] () { return store_results == STORE_RESULTS::USER; }));

	ecfdescription->addSeparator();

	settings = debug = false;
	catalyst = false;
	catalyst_sleep_time = 0;
	REGISTER(settings, ECFMetaData()
			.setdescription({ "Store settings" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(debug, ECFMetaData()
			.setdescription({ "Store FETI related data" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst, ECFMetaData()
			.setdescription({ "In-situ visualization by Catalyst" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst_sleep_time, ECFMetaData()
			.setdescription({ "The sleep time between each time steps when catalyst is used" })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	ecfdescription->addSpace();

	collected = true;
	REGISTER(collected, ECFMetaData()
			.setdescription({ "Collected results" })
			.setdatatype({ ECFDataType::BOOL }));

	ecfdescription->addSpace();

	REGISTER(monitoring, ECFMetaData()
			.setdescription({ "List of monitors", "Monitor" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.setpattern({ "1" }), _root);
}



