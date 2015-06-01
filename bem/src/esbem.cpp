
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sys/time.h>
#include "ProgressMonitor.h"
#include "IOHelper.h"
#include "FullMatrix.h"
#include "SurfaceMesh3D.h"
#include "VolumeMesh3D.h"
#include "BESpace.h"
#include "Vector.h"
#include "BEIntegrator.h"
#include "BEIntegratorLaplace.h"
#include "BEBilinearFormLaplace1Layer.h"
#include "BEBilinearFormLaplace2Layer.h"
#include "BEBilinearFormLaplaceHypersingular.h"
#include "BEIntegratorHelmholtz.h"
#include "BEBilinearFormHelmholtz1Layer.h"
#include "BEBilinearFormHelmholtz2Layer.h"
#include "BEBilinearFormHelmholtzHypersingular.h"
#include "IdentityOperator.h"
#include "LaplaceHypersingularOperator.h"
#include "LaplaceSteklovPoincareOperator.h"
#include "Tree.h"
#include "FastBESpace.h"
#include "RepresentationFormulaLaplace.h"
#include "RepresentationFormulaHelmholtz.h"
#include "SparseMatrix.h"
#include "BlockMatrix.h"
#include "BESpaceTime.h"
#include "BEBilinearFormWave1Layer.h"
#include "BEBilinearFormWaveHypersingular.h"
#include "MPIBlockMatrix.h"
#include "LeftPreconditioner.h"
#include "WavePreconditioner.h"
#include "PotentialsWave.h"
#include "PotentialsLaplace.h"
#include "PotentialsHelmholtz.h"
#include "Problem.h"
#include "WaveScatteringProblem.h"
#include "HomogenizationProblem.h"
#include "ACAMatrix.h"
#include "MathFun.h"
#include "HelmholtzRegularizedExteriorNeumannOperator.h"
#include "BEBilinearFormLame1Layer.h"
#include "BEBilinearFormLame2Layer.h"
#include "BEBilinearFormLameHypersingular.h"

#include "esbem.h"

void bem4i::getLameSteklovPoincare(
    double * Sarray,
    esint nNodes,
    const double * nodes,
    esint nElems,
    const esint * elems,
    double nu,
    double E,
    esint orderNear,
    esint orderFar,
    bool verbose
    ) {

  typedef double SC;
  typedef double SCVT;
  typedef esint LO;

  std::vector< SCVT > nodesv;
  nodesv.assign( nodes, nodes + 3 * nNodes );
  std::vector< LO > elemsv;
  elemsv.assign( elems, elems + 3 * nElems );
  SurfaceMesh3D< LO, SC > mesh( nodesv, elemsv );
  if ( verbose ) mesh.printInfo( );

  quadratureType quadType = SauterSchwab;
  int quadNear[ 4 ] = { orderNear, orderNear, orderNear, orderNear };
  int quadFar[] = { orderFar, orderFar };
  if ( verbose ) std::cout << "Sauter-Schwab near-field order " << orderNear <<
      ", Gauss far-field order " << orderFar << std::endl;

  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  SparseMatrix<LO, SC> T12, T13, T23, T;
  mesh.assembleT( T12, T13, T23, T );
  std::vector< SparseMatrix< LO, SC > * > Tv;
  Tv.push_back( &T12 );
  Tv.push_back( &T13 );
  Tv.push_back( &T23 );
  Tv.push_back( &T );

  if ( verbose ) ProgressMonitor::init( "V, Vlap" );
  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame1Layer< LO, SC > formV( &bespace00, quadNear,
      quadType, quadFar, true );
  formV.setE( E );
  formV.setNu( nu );
  formV.assemble( *V, *Vlap );
  if ( verbose ) ProgressMonitor::step( );

  if ( verbose ) ProgressMonitor::init( "Klap" );
  FullMatrix< LO, SC > * Klap = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formKlap( &bespace10, quadNear,
      quadType, quadFar );
  formKlap.assemble( *Klap );
  if ( verbose ) ProgressMonitor::step( );

  if ( verbose ) ProgressMonitor::init( "K" );
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame2Layer< LO, SC > formK( &bespace10, quadNear,
      quadType, quadFar );
  formK.setE( E );
  formK.setNu( nu );
  formK.assemble( *K, *Vlap, *V, *Klap, Tv );
  if ( verbose ) ProgressMonitor::step( );

  delete Klap;

  if ( verbose ) ProgressMonitor::init( "D" );
  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLameHypersingular< LO, SC > formD( &bespace11, quadNear,
      quadType, quadFar );
  formD.setE( E );
  formD.setNu( nu );
  formD.assemble( *D, *Vlap, *V, Tv );
  if ( verbose ) ProgressMonitor::step( );

  delete Vlap;

  IdentityOperator< LO, SC > id( &bespace10 );
  SparseMatrix< LO, SC > M;
  id.assemble( M );

  Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> * Me = M.getEigenSparseMatrix( );

  for ( LO j = 0; j < Me->outerSize( ); ++j ) {
    typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
    it( *Me, j );
    for (; it; ++it ) {
      K->add( it.row( ), it.col( ), 0.5 * it.value( ) );
      K->add( it.row( ) + nElems, it.col( ) + nNodes, 0.5 * it.value( ) );
      K->add( it.row( ) + 2 * nElems, it.col( ) + 2 * nNodes,
          0.5 * it.value( ) );
    }
  }

  FullMatrix< LO, SC > * VinvK = new FullMatrix< LO, SC >( *K );
  if ( verbose ) ProgressMonitor::init( "Applying inverse of V" );
  V->CholeskiSolve( *VinvK );
  if ( verbose ) ProgressMonitor::step( );

  delete V;

  FullMatrix< LO, SC > S( 3 * nNodes, 3 * nNodes, Sarray );
  S.setAll( 0.0 );

  S.multiply( *K, *VinvK, true, false, 1.0, 0.0 );

  delete VinvK;
  delete K;

  S.add( *D, 1.0 );

  delete D;
}
