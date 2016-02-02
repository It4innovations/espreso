
#include "utils.h"



using std::vector;


#pragma once



class DenseMatrixPack
{

public:

  // Default constructor
  DenseMatrixPack();

  // Constructor
  DenseMatrixPack( long maxNMatrices, long preallocSize, int device = 0 );

  // Copy constructor
  DenseMatrixPack( const DenseMatrixPack& orig );

  // Destructor
  ~DenseMatrixPack();

  // deletes current data and resizes the buffer
  void Resize( long maxNMatrices, long preallocSize );

  // prepares arrays containing offsets etc. (necessary to run in cycle before AddDenseMatrix)
  void PreparePack(
    int i,
    int nRows,
    int nCols,
    bool isPacked
  );

  // Adds dense matrix to the pack (if there is enough space)
  void AddDenseMatrix(
    int i,
    double * matrixData
  );

  void AllocateVectors() {
    mic_x_in = ( double * ) malloc( totalCols * sizeof( double ) );
    mic_y_out = ( double * ) malloc( totalRows * sizeof( double ) );
  }

  // Sends matrices to MIC, preallocates data for input/ouptut vectors
  void CopyToMIC();

  // Multiplies input vectors with matrices in pack on cpu
  void DenseMatsVecs(
    char T_for_transpose_N_for_not_transpose
  );

  // Multiplies input vectors with matrices in pack on mic
  void DenseMatsVecsMIC(
    char T_for_transpose_N_for_not_transpose
  );

  // Multiplies input vectors with matrices in pack on mic - start of async. c.
  void DenseMatsVecsMIC_Start(
    char T_for_transpose_N_for_not_transpose
  );

  // Multiplies input vectors with matrices in pack on mic - sync.
  void DenseMatsVecsMIC_Sync( );

  void SetX(
    long vector,
    long position,
    double value
  );

  void GetY(
    long vector,
    std::SEQ_VECTOR <double> & y
  );

  void SetDevice(
    int device
  );

private:

  // MIC number
  int device;

  // maximum number of matrices in pack
  long maxNMatrices;

  // number of matrices in pack
  long nMatrices;

  // size of preallocated data
  long preallocSize;

  // free space
  long freeSpace;

  // array of matrix values
  double * matrices;

  // array with numbers of rows of individual matrices
  int * rows;

  // array with numbers of cols of individual matrices
  int * cols;

  // total number of matrix rows
  long totalRows;

  // total number of matrix cols
  long totalCols;

  // array of offsets to the beginnings of matrices data
  long * offsets;

  // offsets of rows
  long * rowOffsets;

  // offsets of columns
  long * colOffsets;

  // array of lengths of matrix data
  long * lengths;

  // is i-th matrix symmetric and stored in packed format?
  bool * packed;

  // input buffer on MIC
  double * mic_x_in;

  // output buffer on MIC
  double * mic_y_out;

};
