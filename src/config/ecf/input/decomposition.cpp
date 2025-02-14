
#include "decomposition.h"
#include "config/configuration.hpp"

using namespace espreso;

ParMETISConfiguration::ParMETISConfiguration()
{
    refinement = false;
    REGISTER(refinement, ECFMetaData()
            .setdescription({ "Apply refinement." })
            .setdatatype({ ECFDataType::BOOL }));

    tolerance = 1.05;
    REGISTER(tolerance, ECFMetaData()
            .setdescription({ "Imbalance tolerance." })
            .setdatatype({ ECFDataType::FLOAT }));
}

METISConfiguration::METISConfiguration()
{
    ptype = PTYPE::KWAY;
    REGISTER(ptype, ECFMetaData()
            .setdescription({ "Partition method." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("RB").setdescription("Multilevel recursive bisectioning."))
            .addoption(ECFOption().setname("KWAY").setdescription("Multilevel k-way partitioning.")));

    objtype = OBJTYPE::VOL;
    REGISTER(objtype, ECFMetaData()
            .setdescription({ "Partition method." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("CUT").setdescription("Edge-cut minimization."))
            .addoption(ECFOption().setname("VOL").setdescription("Total communication volume minimization."))
            .addoption(ECFOption().setname("NODE").setdescription("Node minimization.")));

    ctype = CTYPE::SHEM;
    REGISTER(ctype, ECFMetaData()
            .setdescription({ "Partition method." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("RM").setdescription("Random matching."))
            .addoption(ECFOption().setname("SHEM").setdescription("Sorted heavy-edge matching.")));

    iptype = IPTYPE::GROW;
    REGISTER(iptype, ECFMetaData()
            .setdescription({ "Partition method." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("GROW").setdescription("Grows a bisection using a greedy strategy."))
            .addoption(ECFOption().setname("RANDOM").setdescription("Computes a bisection at random followed by a refinement."))
            .addoption(ECFOption().setname("EDGE").setdescription("Derives a separator from an edge cut."))
            .addoption(ECFOption().setname("NODE").setdescription("Grow a bisection using a greedy node-based strategy.")));

    rtype = RTYPE::FM;
    REGISTER(rtype, ECFMetaData()
            .setdescription({ "Partition method." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("FM").setdescription("FM-based cut refinement."))
            .addoption(ECFOption().setname("GREEDY").setdescription("Greedy-based cut and volume refinement."))
            .addoption(ECFOption().setname("SEP2SIDED").setdescription("Two-sided node FM refinement."))
            .addoption(ECFOption().setname("SEP1SIDED").setdescription("One-sided node FM refinement.")));

    ncuts = 1;
    REGISTER(ncuts, ECFMetaData()
            .setdescription({ "Specifies the number of different partitionings that it will compute." })
            .setdatatype({ ECFDataType::INTEGER }));

    nseps = 1;
    REGISTER(nseps, ECFMetaData()
            .setdescription({ "Specifies the number of different separators that it will compute at each level of nested dissection." })
            .setdatatype({ ECFDataType::INTEGER }));

    niter = 20;
    REGISTER(niter, ECFMetaData()
            .setdescription({ "Specifies the number of iterations for the refinement algorithms at each stage of the uncoarsening process." })
            .setdatatype({ ECFDataType::INTEGER }));

    seed = 0;
    REGISTER(seed, ECFMetaData()
            .setdescription({ "Specifies the seed for the random number generator." })
            .setdatatype({ ECFDataType::INTEGER }));

    minconn = 0;
    REGISTER(minconn, ECFMetaData()
            .setdescription({ "Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph." })
            .setdatatype({ ECFDataType::INTEGER }));

    no2hp = 0;
    REGISTER(no2hp, ECFMetaData()
            .setdescription({ "Specifies that the coarsening will not perform any 2–hop matchings when the standard matching approach fails to sufficiently coarsen the graph." })
            .setdatatype({ ECFDataType::INTEGER }));

    contig = 1;
    REGISTER(contig, ECFMetaData()
            .setdescription({ "Specifies that the partitioning routines should try to produce partitions that are contiguous." })
            .setdatatype({ ECFDataType::INTEGER }));

    compress = 0;
    REGISTER(compress, ECFMetaData()
            .setdescription({ "Specifies that the graph should be compressed by combining together vertices that have identical adjacency lists." })
            .setdatatype({ ECFDataType::INTEGER }));

    ccorder = 1;
    REGISTER(ccorder, ECFMetaData()
            .setdescription({ "Specifies if the connected components of the graph should first be identified and ordered separately." })
            .setdatatype({ ECFDataType::INTEGER }));

    dbglvl = 0;
    REGISTER(dbglvl, ECFMetaData()
            .setdescription({ "Specifies the amount of progress/debugging information will be printed during the execution of the algorithms." })
            .setdatatype({ ECFDataType::INTEGER }));

    pfactor = 0;
    REGISTER(pfactor, ECFMetaData()
            .setdescription({ "Specifies the minimum degree of the vertices that will be ordered last." })
            .setdatatype({ ECFDataType::INTEGER }));

    ufactor = 100;
    REGISTER(ufactor, ECFMetaData()
            .setdescription({ "Specifies the maximum allowed load imbalance among the partitions (1 + x) / 1000." })
            .setdatatype({ ECFDataType::INTEGER }));
}

PTScotchConfiguration::PTScotchConfiguration()
{

}

ScotchConfiguration::ScotchConfiguration()
{

}

KaHIPConfiguration::KaHIPConfiguration() {

}

DecompositionConfiguration::DecompositionConfiguration()
{
    parallel_decomposer = ParallelDecomposer::PARMETIS;
    REGISTER(parallel_decomposer, ECFMetaData()
            .setdescription({ "Tool that is used for decomoposition. " })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("NONE").setdescription("Skip parallel decomposition."))
            .addoption(ECFOption().setname("METIS").setdescription("METIS library (called by the root process)."))
            .addoption(ECFOption().setname("PARMETIS").setdescription("ParMETIS library."))
            .addoption(ECFOption().setname("PTSCOTCH").setdescription("PT-Scotch library."))
            .addoption(ECFOption().setname("HILBERT_CURVE").setdescription("Sort elements centers according Hilbert space filling curve.")));

    sequential_decomposer = SequentialDecomposer::METIS;
    REGISTER(sequential_decomposer, ECFMetaData()
            .setdescription({ "Tool that is used for decomoposition. " })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("NONE").setdescription("Skip sequential decomposition."))
            .addoption(ECFOption().setname("METIS").setdescription("METIS library."))
            .addoption(ECFOption().setname("Scotch").setdescription("Scotch library."))
            .addoption(ECFOption().setname("KaHIP").setdescription("KaHIP library.")));

    mesh_duplication = 1;
    REGISTER(mesh_duplication, ECFMetaData()
            .setdescription({ "The number of parallel mesh instances (usefull etc. for harmonic solvers)." })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

    domains = 0;
    REGISTER(domains, ECFMetaData()
            .setdescription({ "Number of domains for each cluster (Keep 0 for automatic decomposition)." })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

    ecfdescription->addSpace();

    force_continuity = false;
    REGISTER(force_continuity, ECFMetaData()
            .setdescription({ "Force continuous decomposition." })
            .setdatatype({ ECFDataType::BOOL }));

    separate_materials = false;
    separate_regions = false;
    separate_etypes = false;
    REGISTER(separate_materials, ECFMetaData()
            .setdescription({ "Decomposition respect materials." })
            .setdatatype({ ECFDataType::BOOL }));
    REGISTER(separate_regions, ECFMetaData()
            .setdescription({ "Decomposition respect regions." })
            .setdatatype({ ECFDataType::BOOL }));
    REGISTER(separate_etypes, ECFMetaData()
            .setdescription({ "Decomposition respect elements types." })
            .setdatatype({ ECFDataType::BOOL }));

    REGISTER(parmetis_options, ECFMetaData()
            .setdescription({ "ParMETIS options." }));
    REGISTER(metis_options, ECFMetaData()
            .setdescription({ "ParMETIS options." }));
    REGISTER(ptscotch_options, ECFMetaData()
            .setdescription({ "PT-Scotch options." }));
    REGISTER(scotch_options, ECFMetaData()
            .setdescription({ "PT-Scotch options." }));
    REGISTER(kahip_options, ECFMetaData()
            .setdescription({ "KaHIP options." }));
}




