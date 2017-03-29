#include "DenseMatrixPack.h"

#include "../../../basis/logging/logging.h"

namespace espreso {

    DenseMatrixPack::DenseMatrixPack(const ESPRESOSolver &configuration)
    : configuration(configuration) {

        this->device = 0;

        this->maxNMatrices = 0;
        this->preallocSize = 0;

        // preallocate data structures
        this->matrices = NULL;
        this->rows = NULL;
        this->cols = NULL;
        this->offsets = NULL;
        this->rowOffsets = NULL;
        this->colOffsets = NULL;
        this->lengths =NULL;
        this->packed = NULL;

        this->nMatrices = 0;
        this->freeSpace = 0;
        this->totalRows = 0;
        this->totalCols = 0;

        this->mic_x_in = NULL;
        this->mic_y_out = NULL;

        this->copiedToMIC = false;
        this->MICratio = 1.0;
        this->elapsedTime = new double[1];
    }

    DenseMatrixPack::DenseMatrixPack(
            const ESPRESOSolver &configuration,
            long maxNMatrices,
            long preallocSize,
            int device
            ): configuration(configuration) {

        this->device = device;

        this->maxNMatrices = maxNMatrices;
        this->preallocSize = preallocSize;

        // preallocate data structures
        this->matrices = ( double * ) malloc( preallocSize * sizeof( double ) );
        this->rows = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
        this->cols = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
        this->offsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->rowOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->colOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->lengths = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->packed = ( bool * ) malloc( maxNMatrices * sizeof( bool ) );

        this->nMatrices = 0;
        this->freeSpace = preallocSize;
        this->offsets[ 0 ] = 0;

        this->totalRows = 0;
        this->totalCols = 0;

        this->mic_x_in = NULL;
        this->mic_y_out = NULL;

        this->copiedToMIC = false;
        this->MICratio = 1.0;
        this->elapsedTime = new double[1];
    }

    void DenseMatrixPack::Resize( long maxNMatrices, long preallocSize ) {

        if (this->preallocSize != 0) {
            free( this->matrices );
            free( this->rows );
            free( this->cols );
            free( this->offsets );
            free( this->rowOffsets );
            free( this->colOffsets );
            free( this->lengths );
            free( this->packed );
        }

        this->maxNMatrices = maxNMatrices;
        this->preallocSize = preallocSize;

        // preallocate data structures
        this->matrices = ( double * ) malloc( preallocSize * sizeof( double ) );
        this->rows = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
        this->cols = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
        this->offsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->rowOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->colOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->lengths = ( long * ) malloc( maxNMatrices * sizeof( long ) );
        this->packed = ( bool * ) malloc( maxNMatrices * sizeof( bool ) );

        this->nMatrices = 0;
        this->freeSpace = preallocSize;
        this->offsets[ 0 ] = 0;
        this->rowOffsets[ 0 ] = 0;
        this->colOffsets[ 0 ] = 0;
        this->totalRows = 0;
        this->totalCols = 0;
    }

    DenseMatrixPack::DenseMatrixPack( const DenseMatrixPack& orig ) {

	this->configuration = configuration;
        this->device = orig.device;

        this->maxNMatrices = orig.maxNMatrices;
        this->preallocSize = orig.preallocSize;

        // preallocate data structures
        this->matrices = orig.matrices;
        this->matrices_mic = orig.matrices_mic;
        this->rows = orig.rows;
        this->cols = orig.cols;
        this->offsets = orig.offsets;
        this->rowOffsets = orig.rowOffsets;
        this->colOffsets = orig.colOffsets;
        this->lengths = orig.lengths;
        this->packed = orig.packed;

        this->nMatrices = orig.nMatrices;
        this->freeSpace = orig.freeSpace;
        this->totalRows = orig.totalRows;
        this->totalCols = orig.totalCols;

        this->mic_x_in = orig.mic_x_in;
        this->mic_y_out = orig.mic_y_out;
        this->MICratio = orig.MICratio;
        this->elapsedTime = orig.elapsedTime;

        this->elapsedTime = new double[1];
        this->elapsedTime[0] = orig.elapsedTime[0];
    }

    DenseMatrixPack::~DenseMatrixPack() {
        // free the MIC's memory
        if (this->copiedToMIC) {
#pragma offload_transfer target(mic:device) if(1) \
            nocopy( rows : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( cols : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( offsets : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( rowOffsets : alloc_if( 0 ) free_if(1) ) \
            nocopy( colOffsets : alloc_if( 0 ) free_if(1) ) \
            nocopy( lengths : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( mic_x_in : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( mic_y_out : alloc_if( 0 ) free_if( 1 ) ) \
            nocopy( elapsedTime : alloc_if( 0 ) free_if( 1 ) )
        }

        if (this->matrices != NULL )
            free( this->matrices );
        if (this->rows != NULL )
            free( this->rows );
        if (this->cols != NULL )
            free( this->cols );
        if (this->offsets != NULL )
            free( this->offsets );
        if (this->rowOffsets != NULL )
            free( this->rowOffsets );
        if (this->colOffsets != NULL )
            free( this->colOffsets );
        if (this->lengths != NULL )
            free( this->lengths );
        if (this->packed != NULL )
            free( this->packed );
        if ( this->mic_x_in != NULL )
            free( this->mic_x_in );
        if (this->mic_y_out != NULL )
            free( this->mic_y_out );
        if (this->elapsedTime != NULL) 
            delete [] elapsedTime;

        this->matrices = NULL;
        this->rows = NULL;
        this->cols = NULL;
        this->offsets = NULL;
        this->colOffsets = NULL;
        this->rowOffsets = NULL;
        this->lengths = NULL;
        this->packed = NULL;
        this->mic_x_in = NULL;
        this->mic_y_out = NULL;
    }

    void DenseMatrixPack::SetDevice(
            int device
            ) {
        this->device = device;
    }

    double * DenseMatrixPack::getMatrixPointer(
            eslocal matrix
            ) {
        return this->matrices + this->offsets[matrix];
    }

    void DenseMatrixPack::PreparePack(
            eslocal i,
            eslocal nRows,
            eslocal nCols,
            bool isPacked
            ) {

        long size = 0;
        if (!isPacked) {
            size = nRows * nCols;
        } else {
            // isPacked => is symmetric
            size = ( ( 1.0 + ( double ) nRows ) * ( ( double ) nRows ) / 2.0 );
        }

        if ( i > maxNMatrices ) {
            ESINFO(ERROR) << "Could not add matrix. Maximum number of matrices stored.";
        } else if ( size > freeSpace ) {
            ESINFO(ERROR) << "Could not add matrix. Not enough allocated memory.";
        }

        this->packed[ i ] = isPacked;
        this->rows[ i ] = nRows;
        this->cols[ i ] = nCols;
        this->rowOffsets[ i ] = this->totalRows;
        this->colOffsets[ i ] = this->totalCols;
        this->totalRows += nRows;
        this->totalCols += nCols;

        this->lengths[ i ] = size;
        if ( i > 0 ) {
            this->offsets[ i ] = this->offsets[ i - 1 ] + this->lengths[ i - 1 ];
        }
        freeSpace -= size;
        this->nMatrices++;
    }

    void DenseMatrixPack::AddDenseMatrix(
            eslocal i,
            double * matrixData
            ) {

        if ( i >= nMatrices ) {
            ESINFO(ERROR) << "Could not add matrix. Maximum number of matrices stored.";
            return;
        }
        long size = this->lengths[i];

        if ( !this->packed[i] ) {
            memcpy( this->matrices + this->offsets[ i ],
                    matrixData, sizeof( double ) * size );
        } else {
            for (long j = 0; j < this->cols[ i ]; j++ ) {
                long pos = ((double)( 1.0 +j ))*((double) j)/2.0;
                memcpy( this->matrices + this->offsets[ i ] + pos,
                        matrixData + j * this->rows[ i ], sizeof( double ) * (j+1) );
            }
        }
    }

    void DenseMatrixPack::SetX(
            long vector,
            long position,
            double value
            ) {
        if ( this->mic_x_in != NULL ) {
            this->mic_x_in[position + this->colOffsets[ vector ]] = value;
        } else {
            ESINFO(ERROR) << "Could not set vector element. Vector not yet allocated.";
        }
    }

    void DenseMatrixPack::GetY(
            long vector,
            std::SEQ_VECTOR <double> & y
            ) {
        if ( this->mic_y_out != NULL ) {
            memcpy( &(y[0]), this->mic_y_out + this->rowOffsets[ vector ],
                    this->rows[vector] * sizeof( double ) );
            //std::cout << this->rows[vector] << " "<< y.size() << std::endl;
            //for (int i = 0 ; i < y.size(); i++) {
            //  std::cout << "nastavuju" << i << std::endl;
            //  y[i] = 1.0;
            //}
        } else {
            ESINFO(ERROR) << "Could not copy vector. Vector not yet allocated.";
        }
    }

    void DenseMatrixPack::CopyToMIC( ) {
#ifdef MIC
        long matrixLength = preallocSize - freeSpace;
        double * tmp_pointer;
        // allocate targetptr array on MIC
#pragma offload_transfer target(mic:device) \
        nocopy(tmp_pointer : length( matrixLength) alloc_if(1) free_if(0) targetptr) \
        in(this : alloc_if(1) free_if(0))

        // preallocated input/output buffers
        //mic_x_in = ( double * ) malloc( totalCols * sizeof( double ) );
        //mic_y_out = ( double * ) malloc( totalRows * sizeof( double ) );
        this->matrices_mic = tmp_pointer;
#pragma offload_transfer target(mic:device) if(1) \
        in( matrices : length( matrixLength ) into(matrices_mic) alloc_if( 0 ) free_if( 0 ) targetptr ) \
        in( rows : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( cols : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( rowOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
        in( colOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
        in( lengths : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( packed : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( mic_x_in : length( totalCols ) alloc_if( 1 ) free_if( 0 ) ) \
        in( mic_y_out : length( totalRows ) alloc_if( 1 ) free_if( 0 ) ) \
        in( elapsedTime : length( 1 ) alloc_if( 1 ) free_if( 0 ) ) \
        in( this : length(0) alloc_if( 0 ) free_if( 0 ) )

#endif
        this->copiedToMIC = true;
        if ( !this->loadBalancing ) {
            free(this->matrices);
            this->matrices = NULL;
        }
    }

    void DenseMatrixPack::DenseMatsVecs(
            char T_for_transpose_N_for_not_transpose
            ) {

        double alpha = 1.0;
        double beta  = 0.0;
        eslocal one = 1;

        omp_set_num_threads(24);

#pragma omp parallel for //schedule(dynamic,10)
        for ( long i = 0 ; i < nMatrices; i++ ) {
            //TODO:uncoment to fix
            //        if ( !packed[i] ) {
            //        dgemv(&T_for_transpose_N_for_not_transpose,
            //  				&(rows[i]), &(cols[i]),
            //  				&alpha, matrices + offsets[i], &(rows[i]),
            //  				mic_x_in + colOffsets[i], &one,
            //  				&beta, mic_y_out + rowOffsets[i], &one);
            //        } else {
            //  				cblas_dspmv(CblasColMajor, CblasUpper,
            //  					rows[i], 1.0, matrices + offsets[i], mic_x_in + colOffsets[i],
            //            1, 0.0, mic_y_out + rowOffsets[i], 1);
            //        }
        }


    }

    void DenseMatrixPack::DenseMatsVecsCPU(
            long start,
            long end,
            char T_for_transpose_N_for_not_transpose
            ) {
        double alpha = 1.0;
        double beta  = 0.0;
        eslocal one = 1;

#pragma omp parallel for //schedule(dynamic,10)
        for ( long i = start ; i < end; i++ ) {
            if ( !packed[i] ) {
                dgemv(&T_for_transpose_N_for_not_transpose,
                        &(rows[i]), &(cols[i]),
                        &alpha, matrices + offsets[i], &(rows[i]),
                        mic_x_in + colOffsets[i], &one,
                        &beta, mic_y_out + rowOffsets[i], &one);
            } else {
                cblas_dspmv(CblasColMajor, CblasUpper,
                        rows[i], 1.0, matrices + offsets[i], mic_x_in + colOffsets[i],
                        1, 0.0, mic_y_out + rowOffsets[i], 1);
            }
        }
    }


    void DenseMatrixPack::DenseMatsVecsRestCPU(
            char T_for_transpose_N_for_not_transpose
            ) {
        double alpha = 1.0;
        double beta  = 0.0;
        eslocal one = 1;
        long start = (long) (MICratio * nMatrices); 
#pragma omp parallel for //schedule(dynamic,10)
        for ( long i = start ; i < nMatrices; i++ ) {
            if ( !packed[i] ) {
                dgemv(&T_for_transpose_N_for_not_transpose,
                        &(rows[i]), &(cols[i]),
                        &alpha, matrices + offsets[i], &(rows[i]),
                        mic_x_in + colOffsets[i], &one,
                        &beta, mic_y_out + rowOffsets[i], &one);
            } else {
                cblas_dspmv(CblasColMajor, CblasUpper,
                        rows[i], 1.0, matrices + offsets[i], mic_x_in + colOffsets[i],
                        1, 0.0, mic_y_out + rowOffsets[i], 1);
            }
        }
    }

    void DenseMatrixPack::DenseMatsVecsMIC(
            char T_for_transpose_N_for_not_transpose
            ) {
#ifdef MIC


#pragma offload target(mic:device) if(1) \
        in( matrices_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
        in( mic_x_in : length( totalCols ) alloc_if( 0 ) free_if( 0 ) ) \
        in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
        in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
        in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        out( mic_y_out : length( totalRows ) alloc_if( 0 ) free_if( 0 ) ) \
        in( this : length(0) alloc_if( 0 ) free_if( 0 ) )
        {
            double alpha = 1.0;
            double beta  = 0.0;
            eslocal one = 1;

            //  omp_set_num_threads(24);

#pragma omp parallel for //schedule(dynamic,10)
            for ( eslocal i = 0 ; i < nMatrices; i++ ) {

                if ( !packed[i] ) {
                    dgemv(&T_for_transpose_N_for_not_transpose,
                            &(rows[i]), &(cols[i]),
                            &alpha, matrices_mic + offsets[i], &(rows[i]),
                            mic_x_in + colOffsets[i], &one,
                            &beta, mic_y_out + rowOffsets[i], &one);
                } else {
                    cblas_dspmv(CblasColMajor, CblasUpper,
                            rows[i], 1.0, matrices_mic + offsets[i], mic_x_in + colOffsets[i],
                            1, 0.0, mic_y_out + rowOffsets[i], 1);
                }
            }
        }

#endif

    }

    void DenseMatrixPack::DenseMatsVecsMIC_Start(
            char T_for_transpose_N_for_not_transpose
            ) {
        long nMatrices = this->nMatrices;
#pragma offload target(mic:device) if(1) \
        in( this->mic_x_in :length(totalCols) alloc_if(0) free_if(0) ) \
        in( matrices_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr) \
        in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
        in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
        in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
        in( mic_y_out : length( 0 ) alloc_if( 0 ) free_if( 0 )  ) \
        in( MICratio ) \
        in( elapsedTime : length(0) alloc_if(0) free_if(0) ) \ 
        signal( mic_y_out )\
            in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
            {
                double alpha = 1.0;
                double beta  = 0.0;
                eslocal one = 1;
                long nIters = (long) (nMatrices*MICratio);
                double start = omp_get_wtime();
                int nth;
#pragma omp parallel for schedule(dynamic)
                for ( long i = 0 ; i < nIters; i++ ) {
                    if ( !packed[i] ) {
                        dgemv(&T_for_transpose_N_for_not_transpose,
                                &(rows[i]), &(cols[i]),
                                &alpha, matrices_mic + offsets[i], &(rows[i]),
                                mic_x_in + colOffsets[i], &one,
                                &beta, mic_y_out + rowOffsets[i], &one);
                    } else {
                        cblas_dspmv(CblasColMajor, CblasUpper,
                                rows[i], 1.0, matrices_mic + offsets[i], mic_x_in + colOffsets[i],
                                1, 0.0, mic_y_out + rowOffsets[i], 1);
                    }
                }
                elapsedTime[0] = omp_get_wtime() - start;
             }
    }

    void DenseMatrixPack::DenseMatsVecsMIC_Sync( ) {
#pragma offload_wait target(mic:device) wait(mic_y_out)
#pragma offload_transfer target(mic:device) \
        out(mic_y_out : length( totalCols ) alloc_if( 0 ) free_if( 0 ) ) \
        out( elapsedTime : length(1) alloc_if(0) free_if(0) )

    }
}
