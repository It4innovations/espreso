
#include "../../generic/utils.h"



using std::vector;


#pragma once

namespace espreso {

class DenseMatrixPack
{
    friend class SparseSolverMIC;

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
    eslocal i,
    eslocal nRows,
    eslocal nCols,
    bool isPacked
  );

  // Adds dense matrix to the pack (if there is enough space)
  void AddDenseMatrix(
    eslocal i,
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

  // Multiplies input vectors with matrices in pack on mic
  void DenseMatsVecsCPU(
    long start,
    long end,
    char T_for_transpose_N_for_not_transpose
  );

  void DenseMatsVecsRestCPU(
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

  int getDevice() {
    return device;
  }

  bool areDataOnMIC() {
    return this->copiedToMIC;  
  }

  double * getMatrixPointer( 
    eslocal matrix
  );

  long getDataLength(
    eslocal matrix
  ) {
    return this->lengths[matrix];
  }

  void setMICratio( 
    double MICratio
   ) {
    this->MICratio = MICratio; 
  }

  double getMICratio() {
    return MICratio;    
  }

long getNMatrices() {
    return this->nMatrices;
  }

double getElapsedTime() {
    return this->elapsedTime;
}

/*
  void setDevice(int device) {
    this->device = device;
  }

  long getMaxNMatrices() {
    return this->maxNMatrices;
  }

  void setMaxNMatrices(long maxNMatrices) {
    this->maxNMatrices = maxNMatrices;
  }

  long getNMatrices() {
    return this->nMatrices;
  }

  void setNMatrices(long nMatrices) {
    this->nMatrices = nMatrices;
  }

  void setPreallocSize(long preallocSize) {
    this->preallocSize = preallocSize;
  }

  long getPreallocSize() {
    return this->preallocSize;
  }

  void setFreeSpace(long freeSpace) {
    this->freeSpace = freeSpace;
  }

  long getFreeSpace() {
    return this->freeSpace;
  }

  void setMatrices(double * matrices) {
    this->matrices = matrices;
  }

  double * getMatrices() {
    return this->matrices;
  }

  void setRows(int * rows) {
    this->rows = rows;
  }

  int * getRows() {
    return this->rows;
  }

  void setCols(int * cols) {
    this->cols = cols;
  }

  int * getCols() {
    return this->cols;
  }

  void setTotalRows(long totalRows) {
    this->totalRows = totalRows;
  }

  long getTotalRows() {
    return this->totalRows;
  }

  void setTotalCols(long totalCols) {
    this->totalCols = totalCols;
  }

  long getTotalCols() {
    return this->totalCols;
  }

  void setOffsets(long * offsets) {
      this->offsets = offsets;
  }

  long * getOffsets() {
    return this->offsets;
  }

  void setRowOffsets(long * rowOffsets) {
      this->rowOffsets = rowOffsets;
  }

  long * getRowOffsets() {
    return this->rowOffsets;
  }

  void setColOffsets(long * colOffsets) {
      this->colOffsets = colOffsets;
  }

  long * getColOffsets() {
    return this->colOffsets;
  }

  void setLengths(long * lenghts) {
    this->lengths = lengths;
  }

  long * getLengths() {
    return this->lengths;
  }
  
  void setPacked(bool * packed) {
    this->packed = packed;
  }

  bool* getPacked() {
    return this->packed;
  }

  void setMic_x_in(double * mic_x_in) {
    this->mic_x_in = mic_x_in;
  }

  double * getMic_x_in() {
    return this->mic_x_in;
  }

  void setMic_y_out(double * mic_y_out) {
    this->mic_y_out = mic_y_out;
  }

  double * getMic_y_out() {
    return this->mic_y_out;
  }

*/
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

  // array of matrix values on MIC
  double * matrices_mic;

  // array with numbers of rows of individual matrices
  eslocal * rows;

  // array with numbers of cols of individual matrices
  eslocal * cols;

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

  // are data copied to MIC
  bool copiedToMIC;

  #pragma offload_attribute(push,target(mic))
  // ratio of work during mv multiplication 
  double MICratio;

  // time for one mv
  double elapsedTime;
#pragma offload_attribute(pop)
};

}
