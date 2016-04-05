/*!
 * @file    MPIBlockACAMatrix.h
 * @author  Michal Merta, Jan Zapletal 
 * @date    February 21, 2016
 * @brief   Header file for the MPIBlockACAMatrix class
 * 
 */

#ifndef MPIBLOCKACAMATRIX_H
#define	MPIBLOCKACAMATRIX_H

#include "MPIBlockMatrix.h"

namespace bem4i {

/*!
 * Class representing block-ACA-sparsified matrix distributed by
 * cyclic graph decomposition
 */
template<class LO, class SC>
class MPIBlockACAMatrix : public MPIBlockMatrix<LO, SC> {
public:

  //! default constructor
  MPIBlockACAMatrix( );

  /*! Constructor taking number of blocks their sizes and rank assignments as arguments
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   * @param[in]   ranks   array of length nBlockRows*nBlockCols containing the assignment of blocks to MPI ranks
   */
  MPIBlockACAMatrix(
      int nBlockRows,
      int nBlockCols,
      LO * numbersOfRows,
      LO * numbersOfCols,
      int * ranks,
      MPI_Comm communicator
      );

  virtual ~MPIBlockACAMatrix( );

  /*!
   * @brief Performs a matrix-vector multiplication
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  virtual void apply(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  void setOuterDOFs( LO i, LO j, std::vector<LO> * dofs ) {
    this->outerDOFs[j * this->nBlockRows + i] = dofs;
  }

  void setInnerDOFs( LO i, LO j, std::vector<LO> * dofs ) {
    this->innerDOFs[j * this->nBlockRows + i] = dofs;
  }
private:



  //! copy constructor
  MPIBlockACAMatrix( const MPIBlockACAMatrix& orig );

  //! 
  std::vector< std::vector<LO>* > outerDOFs;

  std::vector< std::vector<LO>* > innerDOFs;

};

}


#include "MPIBlockACAMatrix.cpp"
#endif	/* MPIBLOCKACAMATRIX_H */

