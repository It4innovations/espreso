
#include "../../generic/utils.h"
#include "../../../configuration/solver/espreso.h"


using std::vector;


#pragma once

namespace espreso {

    class DenseMatrixPack
    {
        friend class SparseSolverMIC;

        public:

        // Default constructor
        DenseMatrixPack( const ESPRESOSolver &configuration);

        // Constructor
        DenseMatrixPack( 
            const ESPRESOSolver &configuration, 
            long maxNMatrices, 
            long preallocSize, 
            int device = 0,
            bool USE_FLOAT = false);

        // Copy constructor
        DenseMatrixPack( const DenseMatrixPack& orig );

        // Destructor
        ~DenseMatrixPack();

        // deletes current data and resizes the buffer
        void Resize( long maxNMatrices, long preallocSize, bool USE_FLOAT = false );

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
            if (!USE_FLOAT) { 
                mic_x_in = ( double * ) malloc( totalCols * sizeof( double ) );
                mic_y_out = ( double * ) malloc( totalRows * sizeof( double ) );
                for (long i = 0 ; i < totalCols; ++i ) {
                    mic_x_in[i] = 0.0;
                }
                for (long i = 0 ; i < totalRows; ++i ) {
                    mic_y_out[i] = 0.0;
                }
            } else {
                mic_x_in_fl = ( float * ) malloc( totalCols * sizeof( float ) );
                mic_y_out_fl = ( float * ) malloc( totalRows * sizeof( float ) );
                for (long i = 0 ; i < totalCols; ++i ) {
                    mic_x_in_fl[i] = 0.0;
                }
                for (long i = 0 ; i < totalRows; ++i ) {
                    mic_y_out_fl[i] = 0.0;
                }

            }

        }

        // Sends matrices to MIC, preallocates data for input/ouptut vectors
        void CopyToMIC();

        // Multiplies input vectors with matrices in pack on mic
        void DenseMatsVecsMIC(
                char T_for_transpose_N_for_not_transpose, 
                double alpha = 1.0, 
                double beta = 0.0
                );

        // Multiplies input vectors with matrices in pack on mic
        void DenseMatsVecsCPU(
                long start,
                long end,
                char T_for_transpose_N_for_not_transpose,
                double alpha = 1.0, 
                double beta = 0.0
                );

        void DenseMatsVecsRestCPU(
                char T_for_transpose_N_for_not_transpose,
                double alpha = 1.0,
                double beta = 0.0
                );

        // Multiplies input vectors with matrices in pack on mic - start of async. c.
        void DenseMatsVecsMIC_Start(
                char T_for_transpose_N_for_not_transpose,
                double alpha = 1.0,
                double beta = 0.0
                );

        // Multiplies input vectors with matrices in pack on mic - sync.
        void DenseMatsVecsMIC_Sync( );

        // sets a value in vector number vector on position position 
        void SetX(
                long vector,
                long position,
                double value
                );

        // returns vector number vector as a reference
        void GetY(
                long vector,
                std::SEQ_VECTOR <double> & y
                );

        // sets a MIC device an instance of this object uses
        void SetDevice(
                int device
                );

        // returns a MIC device an instance of this object uses
        int getDevice() {
            return device;
        }

        // returns true if the necessary data has been copied to MIC
        bool areDataOnMIC() {
            return this->copiedToMIC;  
        }

        // returns a pointer to internal matrix array
        double * getMatrixPointer( 
                eslocal matrix
                );

        // returns a pointer to internal single precision matrix array
        float * getMatrixPointer_fl( 
                eslocal matrix
                );

        // returns lenght of data
        long getDataLength(
                eslocal matrix
                ) {
            return this->lengths[matrix];
        }

        // for load balancing: sets a ratio of work done using MIC
        void setMICratio( 
                double MICratio
                ) {
            this->MICratio = MICratio; 
        }

        // for load balancing: returns a ratio of work done using MIC
        double getMICratio() {
            return MICratio;    
        }

        // returns a number of matrices in the pack
        long getNMatrices() {
            return this->nMatrices;
        }

        // enables load balancing between host and coprocessor
        void enableLoadBalancing() {
            this->loadBalancing = true;
        }

        // disables load balancing between host and coprocessor
        void disableLoadBalancing() {
            this->loadBalancing = false;
        }

        // returns true if load balancing is enabled
        bool getLoadBalancing() {
            return this->loadBalancing;
        }

        // returns time for computation on MIC
        double getElapsedTime() {
            return this->elapsedTime[0];
        }

        // whether to use single precision
        bool USE_FLOAT;

        private:

        ESPRESOSolver configuration;

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

        // array of float matrix values
        float * matrices_fl;

        // array of matrix values on MIC
        double * matrices_mic;

        // array of matrix values on MIC in float
        float * matrices_mic_fl;

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

        // input buffer on MIC
        float * mic_x_in_fl;

        // output buffer on MIC
        float * mic_y_out_fl;


        // are data copied to MIC
        bool copiedToMIC;

        // whether to use load balancing between host and MIC  
        bool loadBalancing;

#pragma offload_attribute(push,target(mic))
        // ratio of work during mv multiplication 
        double MICratio;
        // time for one mv
        double *  elapsedTime;
#pragma offload_attribute(pop)
    };

}
