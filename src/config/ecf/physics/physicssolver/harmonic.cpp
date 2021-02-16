
#include "harmonic.h"
#include "config/configuration.hpp"

espreso::AlternatingFrequencyTime::AlternatingFrequencyTime()
{
	type = TYPE::USER;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("USER").setdescription("User intput."))
			.addoption(ECFOption().setname("AUTOMATIC").setdescription("Automatic.")));

	time_samples = 20;
	REGISTER(time_samples, ECFMetaData()
			.setdescription({ "Number of time samples." })
			.setdatatype({ ECFDataType::INTEGER }));
}


espreso::HarmonicSolverConfiguration::HarmonicSolverConfiguration()
{
	frequency_interval_type = INTERVAL_TYPE::LINEAR;
	REGISTER(frequency_interval_type, ECFMetaData()
			.setdescription({ "Interval type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear spacing."))
			.addoption(ECFOption().setname("LOGARITHMIC").setdescription("Logatirhmic spacing."))
			.addoption(ECFOption().setname("OCTAVE_BAND").setdescription("Octave band.")));

	ecfdescription->addSpace();

	central_frequency = 500;
	REGISTER(central_frequency, ECFMetaData()
			.setdescription({ "Central frequency" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return frequency_interval_type == INTERVAL_TYPE::OCTAVE_BAND; }));

	octave_band = OCTAVE_BAND::BAND_1;
	REGISTER(octave_band, ECFMetaData()
			.setdescription({ "Octave band settings." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("1").setdescription("1"))
			.addoption(ECFOption().setname("1/2").setdescription("1/2"))
			.addoption(ECFOption().setname("1/3").setdescription("1/3"))
			.addoption(ECFOption().setname("1/6").setdescription("1/6"))
			.addoption(ECFOption().setname("1/12").setdescription("1/12"))
			.addoption(ECFOption().setname("1/24").setdescription("1/24"))
			.allowonly([&] () { return frequency_interval_type == INTERVAL_TYPE::OCTAVE_BAND; }));

	min_frequency = 0;
	REGISTER(min_frequency, ECFMetaData()
			.setdescription({ "Central frequency" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return frequency_interval_type != INTERVAL_TYPE::OCTAVE_BAND; }));

	max_frequency = 200;
	REGISTER(max_frequency, ECFMetaData()
			.setdescription({ "Central frequency" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return frequency_interval_type != INTERVAL_TYPE::OCTAVE_BAND; }));

	num_samples = 200;
	REGISTER(num_samples, ECFMetaData()
			.setdescription({ "Central frequency" })
			.setdatatype({ ECFDataType::FLOAT }));
	
	prestress = false;
	REGISTER(prestress, ECFMetaData()
			.setdescription({ "Insert mass matrix of domain interfaces into harmonic matrix." })
			.setdatatype({ ECFDataType::FLOAT }));

	mass_stabilization = false;
	REGISTER(mass_stabilization, ECFMetaData()
			.setdescription({ "Insert mass matrix of domain interfaces into harmonic matrix." })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(damping, ECFMetaData().setdescription({ "Damping configuration." }));
	
	REGISTER(aft, ECFMetaData().setdescription({ "AFT settings." }));

}
