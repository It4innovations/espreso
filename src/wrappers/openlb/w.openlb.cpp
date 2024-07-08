
#include "w.openlb.h"

using namespace espreso;

#ifdef HAVE_OPENLB

#define PARALLEL_MODE_MPI
#define PLATFORM_CPU_SIMD

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

typedef D3Q19<FORCE> NSDESCRIPTOR;
typedef D3Q7<VELOCITY> TDESCRIPTOR;

// Parameters for the simulation setup
T Ra = 1e3;  // Rayleigh-Zahl
const T Pr = 0.71; // Prandtl-Zahl

T lx;

int N = 64; // resolution of the model

const T maxPhysT = 1e4;   // max. simulation time in s, SI unit
const T epsilon = 1.e-3;  // precision of the convergence (residuum)

const T Tcold = 275.15;
const T Thot = 285.15;
const T Tmean = (Tcold + Thot) / 2.0;

/// Values from the literature studies from Davis
T LitVelocity3[] = { 3.649, 3.696, 1.013 };
T LitPosition3[] = { 0.813, 0.178 };
T LitVelocity4[] = { 16.178, 19.617, 1.212 };
T LitPosition4[] = { 0.823, 0.119 };
T LitVelocity5[] = { 34.730, 68.590, 1.975 };
T LitPosition5[] = { 0.855, 0.066 };
T LitVelocity6[] = { 64.530, 219.36, 3.400 };
T LitPosition6[] = { 0.850, 0.036 };
T LitNusselt3 = 1.117;
T LitNusselt4 = 2.238;
T LitNusselt5 = 4.509;
T LitNusselt6 = 8.817;

/// Compute the nusselt number at the left wall
static T computeNusselt(SuperGeometry<T,3>& superGeometry,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice)
{
  int voxel = 0, material = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0;
  T q = 0;

  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); iC++) {
    int ny = NSlattice.getBlock(iC).getNy();

    int iX = 0;
    int iZ = 1;

    for (int iY = 0; iY < ny; ++iY) {
      material = superGeometry.getBlockGeometry(iC).getMaterial(iX,iY,iZ);

      T_x = ADlattice.getBlock(iC).get(iX,iY,iZ).computeRho();
      T_xplus1 = ADlattice.getBlock(iC).get(iX+1,iY,iZ).computeRho();
      T_xplus2 = ADlattice.getBlock(iC).get(iX+2,iY,iZ).computeRho();

      if ( material == 2 ) {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(q, MPI_SUM);
  singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
#endif

  return q / (T)voxel;
}


/// Stores geometry information in form of material numbers
static void prepareGeometry(SuperGeometry<T,3>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,4);

  std::vector<T> extend(3,T());
  extend[0] = lx;
  extend[1] = lx;
  extend[2] = 3.0 * converter.getPhysLength(1);
  std::vector<T> origin(3,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  origin[2] = 0.0;
  IndicatorCuboid3D<T> cuboid2(extend, origin);

  superGeometry.rename(4,1,cuboid2);

  std::vector<T> extendwallleft(3,T(0));
  extendwallleft[0] = converter.getPhysLength(1);
  extendwallleft[1] = lx;
  extendwallleft[2] = 0.1;
  std::vector<T> originwallleft(3,T(0));
  originwallleft[0] = 0.0;
  originwallleft[1] = 0.0;
  originwallleft[2] = 0.0;
  IndicatorCuboid3D<T> wallleft(extendwallleft, originwallleft);

  std::vector<T> extendwallright(3,T(0));
  extendwallright[0] = converter.getPhysLength(1);
  extendwallright[1] = lx;
  extendwallright[2] = 0.1;
  std::vector<T> originwallright(3,T(0));
  originwallright[0] = lx+converter.getPhysLength(1);
  originwallright[1] = 0.0;
  originwallright[2] = 0.0;
  IndicatorCuboid3D<T> wallright(extendwallright, originwallright);

  superGeometry.rename(4,2,1,wallleft);
  superGeometry.rename(4,3,1,wallright);


  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

static void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T omega  = converter.getLatticeRelaxationFrequency();
  T Tomega  = converter.getLatticeThermalRelaxationFrequency();

  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  setBounceBackBoundary(ADlattice, superGeometry, 4);

  NSlattice.defineDynamics<ForcedBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  setBounceBackBoundary(NSlattice, superGeometry, 4);

  /// sets boundary
  setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, superGeometry.getMaterialIndicator({2, 3}));
  setLocalVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry.getMaterialIndicator({2, 3}));

  /// define initial conditions
  AnalyticalConst3D<T,T> rho(1.);
  AnalyticalConst3D<T,T> u0(0.0, 0.0, 0.0);
  AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst3D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst3D<T,T> T_mean(converter.getLatticeTemperature(Tmean));

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry, 1, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 1, rho, u0);
  NSlattice.defineRhoU(superGeometry, 2, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 2, rho, u0);
  NSlattice.defineRhoU(superGeometry, 3, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 3, rho, u0);

  ADlattice.defineRho(superGeometry, 1, T_mean);
  ADlattice.iniEquilibrium(superGeometry, 1, T_mean, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);

  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

static void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                       SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry<T,3>& superGeometry)
{

  // nothing to do here

}

static void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry<T,3>& superGeometry,
                util::Timer<T>& timer,
                bool converged)
{
  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter3D<T> vtkWriter("squareCavity3d");
  SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature3D<T, NSDESCRIPTOR,TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( geometry );
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  AnalyticalFfromSuperF3D<T> interpolation(velocity, true);

  const int vtkIter = 2000.;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT%vtkIter == 0 || converged) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    timer.update(iT);
    timer.printStep();

    /// NSLattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADLattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    vtkWriter.write(iT);
    const double a[3] = {0, 0, 1.};
    BlockReduction3D2D<T> planeReduction(temperature, a);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold*0.98, Thot*1.02, iT, "temperature");

    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction2(normVel, {0, 0, 1});
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

  }

  if ( converged ) {

    T nusselt = computeNusselt(superGeometry, NSlattice, ADlattice);

    /// Initialize vectors for data output
    T xVelocity[3] = { T() };
    T outputVelX[3] = { T() };
    T yVelocity[3] = { T() };
    T outputVelY[3] = { T() };
    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[3] = { lx / 2, lx * n / (T) outputSize, lx / N * 2 / 2 };
      T xPosition[3] = { lx * n / (T) outputSize, lx / 2, lx / N * 2 / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    // compare to De Vahl Davis' benchmark solutions
    clout << "Comparison against De Vahl Davis (1983):" << std::endl;
    if (Ra == 1e3) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity3[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity3[0]) << std::endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity3[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity3[1]) << std::endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) util::fabs((LitVelocity3[2] - outputVelY[0] / outputVelX[0])  / LitVelocity3[2]) << std::endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition3[0] - outputVelX[1] / lx) / LitPosition3[0]) << std::endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition3[1] - outputVelY[1] / lx) / LitPosition3[1]) << std::endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((LitNusselt3 - nusselt) / nusselt) << std::endl;
    }
    else if (Ra == 1e4) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity4[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity4[0]) << std::endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity4[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity4[1]) << std::endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) util::fabs((LitVelocity4[2] - outputVelY[0] / outputVelX[0])  / LitVelocity4[2]) << std::endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition4[0] - outputVelX[1] / lx) / LitPosition4[0]) << std::endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition4[1] - outputVelY[1] / lx) / LitPosition4[1]) << std::endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((LitNusselt4 - nusselt) / nusselt) << std::endl;
    }
    else if (Ra == 1e5) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity5[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity5[0]) << std::endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity5[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity5[1]) << std::endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) util::fabs((LitVelocity5[2] - outputVelY[0] / outputVelX[0])  / LitVelocity5[2]) << std::endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition5[0] - outputVelX[1] / lx) / LitPosition5[0]) << std::endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition5[1] - outputVelY[1] / lx) / LitPosition5[1]) << std::endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((LitNusselt5 - nusselt) / nusselt) << std::endl;
    }
    else if (Ra == 1e6) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity6[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity6[0]) << std::endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) util::fabs((LitVelocity6[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity6[1]) << std::endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) util::fabs((LitVelocity6[2] - outputVelY[0] / outputVelX[0])  / LitVelocity6[2]) << std::endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition6[0] - outputVelX[1] / lx) / LitPosition6[0]) << std::endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) util::fabs((LitPosition6[1] - outputVelY[1] / lx) / LitPosition6[1]) << std::endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((LitNusselt6 - nusselt) / nusselt) << std::endl;
    }
  }
}


bool OpenLB::isLinked()
{
    return true;
}

OpenLB::OpenLB()
: external(nullptr)
{

}

OpenLB::~OpenLB()
{

}

bool OpenLB::set()
{

    return true;
}

bool OpenLB::update()
{

    return true;
}

bool OpenLB::solve()
{
    /// === 1st Step: Initialization ===
    OstreamManager clout(std::cout,"main");
    olbInit(nullptr, nullptr);
    singleton::directories().setOutputDir("./tmp/");

    T tau = 0.9;

//    if (argc>=2) {
//      Ra = atof(argv[1]);
//    }

    lx  = util::pow(Ra * 15.126e-6 * 15.126e-6 / Pr / 9.81 / (Thot - Tcold) / 0.00341, (T) 1/3);  // length of the square
    T charU = 1.0 / lx /( Pr * 25.684e-3 / 15.126e-6 / 1.0 * 1.0 / 25.684e-3);

    if (Ra==1e3) {
      charU *= LitVelocity3[1];
      N = 64;
    }
    if (Ra==1e4) {
      charU *= LitVelocity4[1];
      N = 128;
    }
    if (Ra==1e5) {
      charU *= LitVelocity5[1];
      N = 256;
    }
    if (Ra==1e6) {
      charU *= LitVelocity6[1];
      N = 512;
    }

    ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
      (T) lx / N,
      (T) (tau - 0.5) / descriptors::invCs2<T,NSDESCRIPTOR>() * util::pow((lx/N),2) / 15.126e-6,
      (T) lx,
      (T) charU,
      (T) 15.126e-6,
      (T) 1.0,
      (T) 25.684e-3,
      (T) Pr * 25.684e-3 / 15.126e-6 / 1.0,
      (T) 0.00341,
      (T) Tcold,
      (T) Thot
    );
    converter.print();

    /// === 2nd Step: Prepare Geometry ===
    std::vector<T> extend(3,T());
    extend[0] = lx + 2*converter.getPhysLength(1);
    extend[1] = lx + converter.getPhysLength(1);
    extend[2] = 2.*converter.getPhysLength(1);
    std::vector<T> origin(3,T());
    IndicatorCuboid3D<T> cuboid(extend, origin);

    /// Instantiation of an empty cuboidGeometry
    CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
    cuboidGeometry.setPeriodicity(false, false, true);  // x, y, z

    /// Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    /// Instantiation of a superGeometry
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer);

    prepareGeometry(superGeometry, converter);

    /// === 3rd Step: Prepare Lattice ===

    SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
    SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);

    //prepareLattice and setBoundaryCondition
    prepareLattice(converter, NSlattice, ADlattice, superGeometry);

    T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                                 converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();
    SuperLatticeCoupling coupling(
      NavierStokesAdvectionDiffusionCoupling{},
      names::NavierStokes{}, NSlattice,
      names::Temperature{},  ADlattice);
    coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
      converter.getLatticeTemperature(Tcold));
    coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
      boussinesqForcePrefactor * Vector<T,3>{0.0,1.0,0.0});

    /// === 4th Step: Main Loop with Timer ===
    util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
    timer.start();

    util::ValueTracer<T> converge(6,epsilon);
    for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

      if (converge.hasConverged()) {
        clout << "Simulation converged." << std::endl;
        clout << "Time " << iT << "." << std::endl;

        getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());

        break;
      }

      /// === 5th Step: Definition of Initial and Boundary Conditions ===
      setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

      /// === 6th Step: Collide and Stream Execution ===
      ADlattice.collideAndStream();
      NSlattice.collideAndStream();

      coupling.execute();

      /// === 7th Step: Computation and Output of the Results ===
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
      if (iT % 1000 == 0) {
        converge.takeValue(computeNusselt(superGeometry, NSlattice, ADlattice),true);
      }
    }

    timer.stop();
    timer.printSummary();

    return true;
}

#else

bool OpenLB::isLinked()
{
    return false;
}

OpenLB::OpenLB()
: external(nullptr)
{

}

OpenLB::~OpenLB()
{

}

bool OpenLB::set()
{
    return false;
}

bool OpenLB::update()
{
    return false;
}

bool OpenLB::solve()
{
    return false;
}

#endif
