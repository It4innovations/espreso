#include "DenseMatrixPack.h"

namespace espreso {

    DenseMatrixPack::DenseMatrixPack(const FETISolverConfiguration &configuration)
        : configuration(configuration) {

            this->device = 0;

            this->maxNMatrices = 0;
            this->preallocSize = 0;

            // preallocate data structures
            this->matrices = NULL;
            this->matrices_fl = NULL;
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
            this->mic_x_in_fl = NULL;
            this->mic_y_out_fl = NULL;

            this->copiedToMIC = false;
            this->MICratio = 1.0;
            this->elapsedTime = new double[1];
            this->USE_FLOAT = false;
        }

    DenseMatrixPack::DenseMatrixPack(
            const FETISolverConfiguration &configuration,
            long maxNMatrices,
            long preallocSize,
            int device, 
            bool USE_FLOAT
            ): configuration(configuration) {

        this->device = device;

        this->maxNMatrices = maxNMatrices;
        this->preallocSize = preallocSize;

        // preallocate data structures
        if (!USE_FLOAT) {
            this->matrices = ( double * ) malloc( preallocSize * sizeof( double ) );
        } else {
            this->matrices_fl = ( float * ) malloc( preallocSize * sizeof( float ) );
        }
        this->rows = ( esint * ) malloc( maxNMatrices * sizeof( esint ) );
        this->cols = ( esint * ) malloc( maxNMatrices * sizeof( esint ) );
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
        this->mic_x_in_fl = NULL;
        this->mic_y_out_fl = NULL;

        this->copiedToMIC = false;
        this->MICratio = 1.0;
        this->elapsedTime = new double[1];
        this->USE_FLOAT = USE_FLOAT;
    }

    void DenseMatrixPack::Resize( long maxNMatrices, long preallocSize, bool USE_FLOAT ) {

        this->USE_FLOAT = USE_FLOAT;

        if (this->preallocSize != 0) {
            if (!USE_FLOAT) {
                free( this->matrices );
            } else {
                free( this->matrices_fl );
            }
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
        if ( !USE_FLOAT ) {
            this->matrices = ( double * ) malloc( preallocSize * sizeof( double ) ); 
        } else {
            this->matrices_fl = ( float * ) malloc( preallocSize * sizeof( float ) ); 
        }
        this->rows = ( esint * ) malloc( maxNMatrices * sizeof( esint ) );
        this->cols = ( esint * ) malloc( maxNMatrices * sizeof( esint ) );
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
        this->matrices_fl = orig.matrices_fl;
        this->matrices_mic_fl = orig.matrices_mic_fl;
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
        this->mic_x_in_fl = orig.mic_x_in_fl;
        this->mic_y_out_fl = orig.mic_y_out_fl;
        this->MICratio = orig.MICratio;

        this->elapsedTime = new double[1];
        this->elapsedTime[0] = orig.elapsedTime[0];
        this->USE_FLOAT = orig.USE_FLOAT;
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
            nocopy( elapsedTime : alloc_if( 0 ) free_if( 1 ) )
            if (!USE_FLOAT) {
#pragma offload_transfer target(mic:device) \
                nocopy( mic_x_in : alloc_if( 0 ) free_if( 1 ) ) \
                nocopy( mic_y_out : alloc_if( 0 ) free_if( 1 ) )
            } else {
#pragma offload_transfer target(mic:device) \
                nocopy( mic_x_in_fl : alloc_if( 0 ) free_if( 1 ) ) \
                nocopy( mic_y_out_fl : alloc_if( 0 ) free_if( 1 ) ) 
            }


        }


        if (this->matrices != NULL )
            free( this->matrices );
        if (this->matrices_fl != NULL ) 
            free( this->matrices_fl );
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
        if ( this->mic_x_in_fl != NULL )
            free( this->mic_x_in_fl);
        if (this->mic_y_out_fl != NULL )
            free( this->mic_y_out_fl );

        if (this->elapsedTime != NULL) 
            delete [] elapsedTime;

        this->matrices = NULL;
        this->matrices_fl = NULL;
        this->rows = NULL;
        this->cols = NULL;
        this->offsets = NULL;
        this->colOffsets = NULL;
        this->rowOffsets = NULL;
        this->lengths = NULL;
        this->packed = NULL;
        this->mic_x_in = NULL;
        this->mic_y_out = NULL;
        this->mic_x_in_fl = NULL;
        this->mic_y_out_fl = NULL;

    }

    void DenseMatrixPack::SetDevice(
            int device
            ) {
        this->device = device;
    }

    double * DenseMatrixPack::getMatrixPointer(
            esint matrix
            ) {
        return this->matrices + this->offsets[matrix];
    }

    float * DenseMatrixPack::getMatrixPointer_fl(
            esint matrix
            ) {
        return this->matrices_fl + this->offsets[matrix];
    }

    void DenseMatrixPack::PreparePack(
            esint i,
            esint nRows,
            esint nCols,
            bool isPacked
            ) {

        long size = 0;
        if (!isPacked) {
            size = nRows * nCols;
        } else {
            // isPacked => is symmetric
            size = ( ( 1.0 + ( double ) nRows ) * ( ( double ) nRows ) / 2.0 );
        }

//        if ( i > maxNMatrices ) {
//            ESINFO(ERROR) << "Could not add matrix. Maximum number of matrices stored.";
//        } else if ( size > freeSpace ) {
//            ESINFO(ERROR) << "Could not add matrix. Not enough allocated memory.";
//        }

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
            esint i,
            double * matrixData
            ) {

        if ( i >= nMatrices ) {
//            ESINFO(ERROR) << "Could not add matrix. Maximum number of matrices stored.";
            return;
        }
        long size = this->lengths[i];


        if (!USE_FLOAT) {
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
        } else {
            if ( !this->packed[i] ) {
                memcpy( this->matrices_fl + this->offsets[ i ],
                        matrixData, sizeof( float ) * size );
            } else {
                for (long j = 0; j < this->cols[ i ]; j++ ) {
                    long pos = ((double)( 1.0 +j ))*((double) j)/2.0;
                    memcpy( this->matrices_fl + this->offsets[ i ] + pos,
                            matrixData + j * this->rows[ i ], sizeof( float ) * (j+1) );
                }
            }

        }
    }

    void DenseMatrixPack::SetX(
            long vector,
            long position,
            double value
            ) {
        if ( !USE_FLOAT ) {
            this->mic_x_in[position + this->colOffsets[ vector ]] = value;
        } else {
            this->mic_x_in_fl[position + this->colOffsets[ vector ]] = (float) value;
        }
    }

    void DenseMatrixPack::GetY(
            long vector,
            std::SEQ_VECTOR <double> & y
            ) {
        if ( this->mic_y_out != NULL || this->mic_y_out_fl != NULL) {
            if (!USE_FLOAT) {
                memcpy( &(y[0]), this->mic_y_out + this->rowOffsets[ vector ],
                        this->rows[vector] * sizeof( double ) );
            } else {
                for (esint i = 0 ; i < this->rows[vector]; ++i) {
                    y[i] = (double) this->mic_y_out_fl[i + this->rowOffsets[ vector ]];
                }
            }
        } else {
//            ESINFO(ERROR) << "Could not copy vector. Vector not yet allocated.";
        }
    }

    void DenseMatrixPack::CopyToMIC( ) {

        long matrixLength = preallocSize - freeSpace;
        if (!USE_FLOAT) {
            double * tmp_pointer;
            // allocate targetptr array on MIC
#pragma offload_transfer target(mic:device) \
            nocopy(tmp_pointer : length( matrixLength) alloc_if(1) free_if(0) targetptr) \
            in(this : alloc_if(1) free_if(0))

            this->matrices_mic = tmp_pointer;

#pragma offload_transfer target(mic:device) \
            in( matrices : length( matrixLength ) into(matrices_mic) alloc_if( 0 ) free_if( 0 ) targetptr ) \
            in( mic_x_in : length( totalCols ) alloc_if( 1 ) free_if( 0 ) ) \
            in( mic_y_out : length( totalRows ) alloc_if( 1 ) free_if( 0 ) )
        } else {
            float * tmp_pointer;
            // allocate targetptr array on MIC
#pragma offload_transfer target(mic:device) \
            nocopy(tmp_pointer : length( matrixLength) alloc_if(1) free_if(0) targetptr) \
            in(this : alloc_if(1) free_if(0))

            this->matrices_mic_fl = tmp_pointer;

#pragma offload_transfer target(mic:device) \
            in( matrices_fl : length( matrixLength ) into(matrices_mic_fl) alloc_if( 0 ) free_if( 0 ) targetptr ) \
            in( mic_x_in_fl : length( totalCols ) alloc_if( 1 ) free_if( 0 ) ) \
            in( mic_y_out_fl : length( totalRows ) alloc_if( 1 ) free_if( 0 ) )
        }

#pragma offload_transfer target(mic:device) if(1) \
        in( rows : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( cols : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( rowOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
        in( colOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
        in( lengths : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( packed : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
        in( elapsedTime : length( 1 ) alloc_if( 1 ) free_if( 0 ) ) \
        in( this : length(0) alloc_if( 0 ) free_if( 0 ) )

        // in the case of load balancing the cpu arrays have to be kept
        this->copiedToMIC = true;
        if ( !this->loadBalancing ) {
            if (!USE_FLOAT) {
                free(this->matrices);
                this->matrices = NULL;
            } else {
                free(this->matrices_fl);
                this->matrices_fl = NULL;
            }
        }
    }

    void DenseMatrixPack::DenseMatsVecsCPU(
            long start,
            long end,
            char T_for_transpose_N_for_not_transpose,
            double alpha,
            double beta
            ) {
        esint one = 1;
        CBLAS_TRANSPOSE trans = CblasNoTrans;
        if (T_for_transpose_N_for_not_transpose == 'T') {
            trans = CblasTrans;
        }

#pragma omp parallel for 
        for ( long i = start ; i < end; i++ ) {
            if ( !packed[i] ) {
                if (!USE_FLOAT) {
                    cblas_dgemv(CblasColMajor, trans, rows[i], cols[i], 
                            alpha, matrices + offsets[i], rows[i],
                            mic_x_in + colOffsets[i], 1, 
                            beta, mic_y_out + rowOffsets[i], 1);
                    //                    dgemv(&T_for_transpose_N_for_not_transpose,
                    //                            &(rows[i]), &(cols[i]),
                    //                            &alpha, matrices + offsets[i], &(rows[i]),
                    //                            mic_x_in + colOffsets[i], &one,
                    //                            &beta, mic_y_out + rowOffsets[i], &one);
                } else {
                    cblas_sgemv(CblasColMajor, trans, rows[i], cols[i], 
                            (float) alpha, matrices_fl + offsets[i], rows[i],
                            mic_x_in_fl + colOffsets[i], 1, 
                            (float) beta, mic_y_out_fl + rowOffsets[i], 1);
                }
            } else {
                if (!USE_FLOAT) {
                    cblas_dspmv(CblasColMajor, CblasUpper,
                            rows[i], (float) alpha, matrices + offsets[i], mic_x_in + colOffsets[i],
                            1, (float) beta, mic_y_out + rowOffsets[i], 1);
                } else {
                    cblas_sspmv(CblasColMajor, CblasUpper,
                            rows[i],(float) alpha, matrices_fl + offsets[i], mic_x_in_fl + colOffsets[i],
                            1, (float) beta, mic_y_out_fl + rowOffsets[i], 1);

                }
            }
        }
    }


    void DenseMatrixPack::DenseMatsVecsRestCPU(
            char T_for_transpose_N_for_not_transpose,
            double alpha,
            double beta
            ) {
        esint one = 1;
        CBLAS_TRANSPOSE trans = CblasNoTrans;
        if (T_for_transpose_N_for_not_transpose == 'T') {
            trans = CblasTrans;
        }

        long start = (long) (MICratio * nMatrices); 
#pragma omp parallel for // schedule(dynamic,1)
        for ( long i = start ; i < nMatrices; i++ ) {
            if ( !packed[i] ) {
                if (!USE_FLOAT) {
                    cblas_dgemv(CblasColMajor, trans, rows[i], cols[i], 
                            alpha, matrices + offsets[i], rows[i],
                            mic_x_in + colOffsets[i], 1, 
                            beta, mic_y_out + rowOffsets[i], 1);

                    //                    dgemv(&T_for_transpose_N_for_not_transpose,
                    //                        &(rows[i]), &(cols[i]),
                    //                        &alpha, matrices + offsets[i], &(rows[i]),
                    //                        mic_x_in + colOffsets[i], &one,
                    //                        &beta, mic_y_out + rowOffsets[i], &one);
                } else {
                    cblas_sgemv(CblasColMajor, trans, rows[i], cols[i], 
                            (float) alpha, matrices_fl + offsets[i], rows[i],
                            mic_x_in_fl + colOffsets[i], 1, 
                            (float) beta, mic_y_out_fl + rowOffsets[i], 1);

                }
            } else {
                if (!USE_FLOAT) {
                    cblas_dspmv(CblasColMajor, CblasUpper,
                            rows[i], (float) alpha, matrices + offsets[i], mic_x_in + colOffsets[i],
                            1, (float) beta, mic_y_out + rowOffsets[i], 1);
                } else {
                    cblas_sspmv(CblasColMajor, CblasUpper,
                            rows[i],(float) alpha, matrices_fl + offsets[i], mic_x_in_fl + colOffsets[i],
                            1, (float) beta, mic_y_out_fl + rowOffsets[i], 1);

                }
            }
        }
    }

    void DenseMatrixPack::DenseMatsVecsMIC(
            char T_for_transpose_N_for_not_transpose, 
            double alpha, 
            double beta
            ) {
        long nMatrices = this->nMatrices;
        long inSize = 0;
        long outSize = 0;
        for ( esint i = 0 ; i < (long) (nMatrices*MICratio); ++i ) {
            inSize += cols[i];
            outSize += rows[i];
        }

        CBLAS_TRANSPOSE trans = CblasNoTrans;
        if (T_for_transpose_N_for_not_transpose == 'T') {
            trans = CblasTrans;
        }

        if (!USE_FLOAT) {

            //            if (beta != 0.0) {
#pragma offload_transfer target(mic:device) \
            in( mic_y_out : length(outSize) alloc_if( 0 ) free_if( 0 ) )
            //          }

#pragma offload target(mic:device) \
            in( this->mic_x_in :length(inSize) alloc_if(0) free_if(0) ) \
            in( matrices_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr) \
            in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( MICratio ) \
            in( elapsedTime : length(0) alloc_if(0) free_if(0) ) \
            in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            out(mic_y_out : length( outSize ) alloc_if( 0 ) free_if( 0 ) ) \
            out( elapsedTime : length(1) alloc_if(0) free_if(0) )
            {
                esint one = 1;
                long nIters = (long) (nMatrices*MICratio);
                double start = omp_get_wtime();
                int nth;
#pragma omp parallel for //schedule(dynamic,1)
                for ( long i = 0 ; i < nIters; i++ ) {
                    if ( !packed[i] ) {
                        cblas_dgemv(CblasColMajor, trans, rows[i], cols[i], 
                                alpha, matrices_mic + offsets[i], rows[i],
                                mic_x_in + colOffsets[i], 1, 
                                beta, mic_y_out + rowOffsets[i], 1);

                        //                                            dgemv(&T_for_transpose_N_for_not_transpose,
                        //                                                    &(rows[i]), &(cols[i]),
                        //                                                    &alpha, matrices_mic + offsets[i], &(rows[i]),
                        //                                                    mic_x_in + colOffsets[i], &one,
                        //                                                    &beta, mic_y_out + rowOffsets[i], &one);
                    } else {
                        cblas_dspmv(CblasColMajor, CblasUpper,
                                rows[i], alpha, matrices_mic + offsets[i], mic_x_in + colOffsets[i],
                                1, beta, mic_y_out + rowOffsets[i], 1);
                    }
                }
                elapsedTime[0] = omp_get_wtime() - start;
            }
        } else {

            if (beta != 0.0) {
#pragma offload_transfer target(mic:device) \
                in( mic_y_out_fl : length(outSize) alloc_if( 0 ) free_if( 0 ) )
            }

#pragma offload target(mic:device) \
            in( this->mic_x_in_fl :length(inSize) alloc_if(0) free_if(0) ) \
            in( matrices_mic_fl : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr) \
            in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( MICratio ) \
            in( elapsedTime : length(0) alloc_if(0) free_if(0) ) \
            in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            out(mic_y_out_fl : length( outSize ) alloc_if( 0 ) free_if( 0 ) ) \
            out( elapsedTime : length(1) alloc_if(0) free_if(0) )
            {
                esint one = 1;
                long nIters = (long) (nMatrices*MICratio);
                double start = omp_get_wtime();
                int nth;
#pragma omp parallel for //schedule(dynamic,1)
                for ( long i = 0 ; i < nIters; i++ ) {
                    if ( !packed[i] ) {
                        cblas_sgemv(CblasColMajor, trans, rows[i], cols[i], 
                                (float) alpha, matrices_mic_fl + offsets[i], rows[i],
                                mic_x_in_fl + colOffsets[i], 1, 
                                (float) beta, mic_y_out_fl + rowOffsets[i], 1);

                        //                    dgemv(&T_for_transpose_N_for_not_transpose,
                        //                            &(rows[i]), &(cols[i]),
                        //                            &alpha, matrices_mic + offsets[i], &(rows[i]),
                        //                            mic_x_in + colOffsets[i], &one,
                        //                            &beta, mic_y_out + rowOffsets[i], &one);
                    } else {
                        cblas_sspmv(CblasColMajor, CblasUpper,
                                rows[i], alpha, matrices_mic_fl + offsets[i], mic_x_in_fl + colOffsets[i],
                                1, beta, mic_y_out_fl + rowOffsets[i], 1);
                    }
                }
                elapsedTime[0] = omp_get_wtime() - start;
            }
        }
    }

    void DenseMatrixPack::DenseMatsVecsMIC_Start(
            char T_for_transpose_N_for_not_transpose,
            double alpha,
            double beta
            ) {
        long nMatrices = this->nMatrices;
        long inSize = 0;
        long outSize =  0;
        for ( esint i = 0 ; i < (long) (nMatrices*MICratio); ++i ) {
            inSize += cols[i];
            outSize += rows[i];
        }

        CBLAS_TRANSPOSE trans = CblasNoTrans;
        if (T_for_transpose_N_for_not_transpose == 'T') {
            trans = CblasTrans;
        }



        if (!USE_FLOAT) {
            if ( beta!=0.0 ) {
#pragma offload_transfer target(mic:device) \
                out( mic_y_out : length(outSize) alloc_if( 0 ) free_if( 0 ) )
            }

#pragma offload target(mic:device) signal(this->mic_y_out) \
            in( this->mic_x_in :length(inSize) alloc_if(0) free_if(0) ) \
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
                in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
                {
                    esint one = 1;
                    long nIters = (long) (nMatrices*MICratio);
                    double start = omp_get_wtime();
                    int nth;
#pragma omp parallel for //schedule(dynamic,1)
                    for ( long i = 0 ; i < nIters; i++ ) {
                        if ( !packed[i] ) {
                            cblas_dgemv(CblasColMajor, trans, rows[i], cols[i], 
                                    alpha, matrices_mic + offsets[i], rows[i],
                                    mic_x_in + colOffsets[i], 1, 
                                    beta, mic_y_out + rowOffsets[i], 1);
                        } else {
                            cblas_dspmv(CblasColMajor, CblasUpper,
                                    rows[i], alpha, matrices_mic + offsets[i], mic_x_in + colOffsets[i],
                                    1, beta, mic_y_out + rowOffsets[i], 1);
                        }
                    }
                    elapsedTime[0] = omp_get_wtime() - start;
                }
        } else {
            if ( beta!=0.0 ) {
#pragma offload_transfer target(mic:device) \
                out( mic_y_out_fl : length(outSize) alloc_if( 0 ) free_if( 0 ) )
            }

#pragma offload target(mic:device) signal(this->mic_y_out) \
            in( this->mic_x_in_fl :length(inSize) alloc_if(0) free_if(0) ) \
            in( matrices_mic_fl : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr) \
            in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
            in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
            in( mic_y_out_fl : length( 0 ) alloc_if( 0 ) free_if( 0 )  ) \
            in( MICratio ) \ 
            in( elapsedTime : length(0) alloc_if(0) free_if(0) ) \
                in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
                {
                    esint one = 1;
                    long nIters = (long) (nMatrices*MICratio);
                    double start = omp_get_wtime();
                    int nth;
#pragma omp parallel for //schedule(dynamic,1)
                    for ( long i = 0 ; i < nIters; i++ ) {
                        if ( !packed[i] ) {
                            cblas_sgemv(CblasColMajor, trans, rows[i], cols[i], 
                                    (float) alpha, matrices_mic_fl + offsets[i], rows[i],
                                    mic_x_in_fl + colOffsets[i], 1, 
                                    (float) beta, mic_y_out_fl + rowOffsets[i], 1);

                        } else {
                            cblas_sspmv(CblasColMajor, CblasUpper,
                                    rows[i], alpha, matrices_mic_fl + offsets[i], mic_x_in_fl + colOffsets[i],
                                    1, beta, mic_y_out_fl + rowOffsets[i], 1);
                        }
                    }
                    elapsedTime[0] = omp_get_wtime() - start;
                }



        }
    }

    void DenseMatrixPack::DenseMatsVecsMIC_Sync( ) {

        long inSize = 0;
        long outSize = 0;
        for ( esint i = 0 ; i < (long) (nMatrices*MICratio); ++i ) {
            outSize += cols[i];
        }
        if (!USE_FLOAT) {
#pragma offload_transfer target(mic:device) wait(this->mic_y_out) \
            out(mic_y_out : length( outSize ) alloc_if( 0 ) free_if( 0 ) ) \
            out( elapsedTime : length(1) alloc_if(0) free_if(0) )
        } else {
#pragma offload_transfer target(mic:device) wait(this->mic_y_out_fl) \
            out(mic_y_out_fl : length( outSize ) alloc_if( 0 ) free_if( 0 ) ) \
            out( elapsedTime : length(1) alloc_if(0) free_if(0) )

        }

    }
}
