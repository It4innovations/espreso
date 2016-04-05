/*!
 * @file    MPIBlockACAMatrix.cpp
 * @author  Michal Merta 
 * @date    November 25, 2013
 * 
 */
#ifdef MPIBLOCKACAMATRIX_H

namespace bem4i {

template<class LO, class SC>
MPIBlockACAMatrix<LO, SC>::MPIBlockACAMatrix( ) {
}

template<class LO, class SC>
MPIBlockACAMatrix<LO, SC>::MPIBlockACAMatrix( const MPIBlockACAMatrix& orig ) {
}

template<class LO, class SC>
MPIBlockACAMatrix<LO, SC>::~MPIBlockACAMatrix( ) {
}

template<class LO, class SC>
MPIBlockACAMatrix<LO, SC>::MPIBlockACAMatrix(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols,
    int * ranks,
    MPI_Comm communicator
    ) {

  // allocate memory for pointers
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  this->numbersOfRows = new LO[nBlockRows];
  this->numbersOfCols = new LO[nBlockCols];
  this->delBlocks = new bool[nBlockRows * nBlockCols];
  this->ranks = new int[nBlockRows * nBlockCols];
  this->nBlockRows = nBlockRows;
  this->nBlockCols = nBlockCols;
  this->communicator = communicator;
  memset( this->delBlocks, 0, nBlockRows * nBlockCols * sizeof ( bool ) );

  // copy numbers of rows and columns of each block
  memcpy( this->numbersOfRows, numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, numbersOfCols, nBlockCols * sizeof ( LO ) );

  memcpy( this->ranks, ranks, nBlockRows * nBlockCols * sizeof ( int ) );

  for ( int i = 0; i < this->nBlockRows * this->nBlockCols; i++ ) {
    this->blocks[i] = nullptr;
    this->outerDOFs[i] = nullptr;
    this->innerDOFs[i] = nullptr;
  }

  this->nRows = 0;
  this->nCols = 0;
  for ( int i = 0; i < this->nBlockRows; i++ ) {
    this->nRows += this->numbersOfRows[i];
  }

  for ( int i = 0; i < this->nBlockCols; i++ ) {
    this->nCols += this->numbersOfCols[i];
  }

}

template<class LO, class SC>
void MPIBlockACAMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  MPI_Barrier( this->communicator );

  // iterate through matrix blocks and multiply
  Vector<LO, SC> yi;
  Vector<LO, SC> xj;
  LO yLength = transA ? this->nCols : this->nRows;

  Vector<LO, SC> *localY = new Vector<LO, SC>( yLength );
  Vector<LO, SC> *globalY = new Vector<LO, SC>( yLength );

  localY->setAll( 0.0 );

  LO maxRows = this->getMaxRows( );
  LO maxCols = this->getMaxCols( );
  SC * yiData;
  SC * xjData;

  if ( !transA ) {
    yiData = new LO[maxRows];
    xjData = new LO[maxCols];
    for ( int i = 0; i < this->nBlockRows; i++ ) {
      for ( int j = 0; j < this->nBlockCols; j++ ) {
        if ( this->getBlock( i, j ) && this->amIOwner( i, j ) ) {
          yi.setData( this->numbersOfRows[i], yiData, false );
          yi.setAll( 0.0 );
          xj.setData( this->numbersOfCols[j], xjData, false );
          for ( LO k = 0; k < this->numbersOfCols[j]; ++k ) {
            xj.set( k, x.get( innerDOFs[j * this->nBlockRows + i]->at( k ) ) );
          }
          this->getBlock( i, j )->apply( *xj, *yi, transA, alpha, 1.0 );
          for ( LO k = 0; k < this->numbersOfRows[i]; ++k ) {
            localY.add( outerDOFs[j * this->nBlockRows + i]->at( k ), yi.get( k ) );
          }
        }
      }
    }
  } else {
    yiData = new LO[maxCols];
    xjData = new LO[maxRows];
    for ( int i = 0; i < this->nBlockCols; i++ ) {
      for ( int j = 0; j < this->nBlockRows; j++ ) {
        if ( this->getBlock( j, i ) && this->amIOwner( j, i ) ) {
          yi.setData( this->numbersOfCols[i], yiData, false );
          yi.setAll(0.0);
          xj.setData( this->numbersOfRows[j], xjData, false );
          for ( LO k = 0; k < this->numbersOfRows[j]; ++k ) {
            xj.set( k, x.get( outerDOFs[i * this->nBlockRows + i]->at( k ) ) );
          }
          this->getBlock( j, i )->apply( *xj, *yi, transA, alpha, 1.0 );
          for ( LO k = 0; k < this->numbersOfCols[i]; ++k ) {
            localY.add( innerDOFs[j * this->nBlockRows + i]->at( k ), yi.get( k ) );
          }
        }
      }
    }
  }
  MPI_Barrier( this->communicator );
  MPI_Allreduce( localY->getData( ), globalY->getData( ), yLength,
      GetType<LO, SC>::MPI_SC( ), MPI_SUM, this->communicator );

  y.scale( beta );
  y.add( *globalY );

  delete localY;
  delete globalY;
  delete yiData;
  delete xjData;

}


}

#endif
