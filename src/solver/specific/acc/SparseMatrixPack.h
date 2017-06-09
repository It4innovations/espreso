
#include "../../generic/utils.h"
#include "../../generic/SparseMatrix.h"
#include "mkl_pardiso.h"


using std::vector;


#pragma once

namespace espreso {

    class SparseMatrixPack
    {
        friend class SparseSolverMIC;

        public:

        // Default constructor
        SparseMatrixPack( const ESPRESOSolver &configuration, bool USE_FLOAT = false);

        // Constructor
        SparseMatrixPack( const ESPRESOSolver &configuration, long maxNMatrices, int device = 0 );

        // Copy constructor
        SparseMatrixPack( const SparseMatrixPack& orig );

        // Destructor
        ~SparseMatrixPack();

        // deletes current data and resizes the buffer
        // void Resize( long maxNMatrices, long preallocSize );

        // prepares arrays containing offsets etc. (necessary to run in cycle before AddDenseMatrix)
        /*void PreparePack(
          eslocal i,
          eslocal nRows,
          eslocal nCols
          bool isPacked
          );

        // Adds sparse matrix to the pack (if there is enough space)
        void AddSparseMatrix(
        eslocal *rowInd,
        eslocal *colInd,
        double * matrixData
        );*/

        // adds n matrices from the array A to the matrix pack
        void AddMatrices(
                SparseMatrix ** A,
                eslocal n,
                eslocal device
                );

        void AllocateVectors() {
            if (!USE_FLOAT) {
                mic_x_in = ( double * ) _mm_malloc( x_in_dim * sizeof( double ), 64 );
                mic_y_out = ( double * ) _mm_malloc( y_out_dim * sizeof( double ), 64 );
            } else {
                mic_x_in_fl = ( float * ) _mm_malloc( x_in_dim * sizeof( float ), 64 );
                mic_y_out_fl = ( float * ) _mm_malloc( y_out_dim * sizeof( float ), 64 );
               }
        }

        // Sends matrices to MIC, preallocates data for input/ouptut vectors
        void CopyToMIC();

        // Solves the system using matices in pack on mic - start of async. c.
        void SolveMIC_Start( );

        // Factorizes matrices copied to the coprocessor
        void FactorizeMIC( );

        // Solves the system using matrices in pack on mic - sync.
        void SolveMIC_Sync( );

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

        int GetDevice() {
            return device;
        }

        bool areDataOnMIC() {
            return this->copiedToMIC;  
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

        void enableLoadBalancing() {
            this->loadBalancing = true;
        }
        void disableLoadBalancing() {
            this->loadBalancing = false;
        }
        bool getLoadBalancing() {
            return this->loadBalancing;
        }

        double getElapsedTime() {
            return this->elapsedTime[0];
        }

        private:

#pragma offload_attribute(push,target(mic))

        ESPRESOSolver configuration;

        // MIC number
        int device;

        // maximum number of matrices in pack
        long maxNMatrices;

        // number of matrices in pack
        long nMatrices;

        // total size of array with matrix values
        long preallocSize;

        // array of matrix values
        double * matrix_values;

        // array of matrix values on MIC (targetptr)
        double * matrix_values_mic;

        // array of matrix values in float
        float * matrix_values_fl;

        // array of matrix values on MIC in float (targetptr)
        float * matrix_values_mic_fl;

        // array of indices of rows
        eslocal * rowInd;

        // array of indices of rows on MIC (targetptr)
        eslocal * rowInd_mic;

        // array of indices of columns
        eslocal * colInd;

        // array of indices of columns on MIC (targetptr)
        eslocal * colInd_mic;

        // array with numbers of rows of individual matrices
        eslocal * rows;

        // array with numbers of cols of individual matrices
        eslocal * cols;

        // number of nonzero entries of individual matrices
        eslocal * nnz;

        // total number of matrix rows
        long totalRows;

        // total number of matrix cols
        long totalCols;

        // domain dim
        long x_in_dim;

        // range dim
        long y_out_dim;

        // array of offsets to the beginnings of matrices data
        long * offsets;

        // offsets of rows
        long * rowOffsets;

        // offsets of columns
        long * colOffsets;
        
        // offsets of input vectors
        long * x_in_offsets;

        // offsets of output vectors
        long * y_out_offsets;

        // input buffer on MIC
        double * mic_x_in;

        // output buffer on MIC
        double * mic_y_out;

        // input buffer on MIC in float
        float * mic_x_in_fl;

        // output buffer on MIC in float
        float * mic_y_out_fl;

        // are data copied to MIC
        bool copiedToMIC;

        // whether to use load balancing between host and MIC  
        bool loadBalancing;

        /* Factorization data */

        MKL_INT mtype;
        
        void *** pt;

        MKL_INT ** iparm;

        double ** dparm;

        MKL_INT ** perm;

        MKL_INT * error;

        // ratio of work during mv multiplication 
        double MICratio;
        // time for one mv
        double *  elapsedTime;

        // whether to use float for computation
        bool USE_FLOAT;

#pragma offload_attribute(pop)
    };

}
