#include "SparseMatrixPack.h"

namespace espreso {

    SparseMatrixPack::SparseMatrixPack(const ESPRESOSolver &configuration)
    : configuration(configuration) {

        this->device = 0;

        this->maxNMatrices = 0;
        this->preallocSize = 0;

        // preallocate data structures
        this->matrix_values = NULL;
        this->matrix_values_mic = NULL;
        this->rows = NULL;
        this->cols = NULL;
        this->rowInd = NULL;
        this->rowInd_mic = NULL;
        this->colInd = NULL;
        this->colInd_mic = NULL;
        this->nnz = NULL;
        this->offsets = NULL;
        this->rowOffsets = NULL;
        this->colOffsets = NULL;

        this->nMatrices = 0;
        this->totalRows = 0;
        this->totalCols = 0;
        this->x_in_dim = 0;
        this->y_out_dim = 0;

        this->mic_x_in = NULL;
        this->mic_y_out = NULL;

        this->copiedToMIC = false;
        this->MICratio = 1.0;
        this->elapsedTime = new double[1];
    }

    SparseMatrixPack::SparseMatrixPack(
            const ESPRESOSolver &configuration,
            long maxNMatrices,
            int device
            ): configuration(configuration) {

        this->device = device;

        this->maxNMatrices = maxNMatrices;
        this->preallocSize = 0;

        this->matrix_values = NULL;
        this->matrix_values_mic = NULL;
        this->rowInd = NULL;
        this->rowInd_mic = NULL;
        this->colInd = NULL;
        this->colInd_mic = NULL;

        // preallocate data structures
        this->rows = NULL;
        this->cols = NULL;
        this->offsets = NULL;
        this->rowOffsets = NULL;
        this->colOffsets = NULL;
        this->nnz = NULL;

        this->nMatrices = 0;
        this->offsets[ 0 ] = 0;
        this->rowOffsets[ 0 ] = 0;
        this->colOffsets[ 0 ] = 0; 

        this->totalRows = 0;
        this->totalCols = 0;
        this->x_in_dim = 0;
        this->y_out_dim = 0;

        this->mic_x_in = NULL;
        this->mic_y_out = NULL;

        this->copiedToMIC = false;
        this->MICratio = 1.0;
        this->elapsedTime = new double[1];
    }

    /*
       void SparseMatrixPack::Resize( long maxNMatrices ) {

       if (this->preallocSize != 0) {
       free( this->rows );
       free( this->cols );
       free( this->offsets );
       free( this->rowOffsets );
       free( this->colOffsets );
       free( this->nnz );
       }

       this->maxNMatrices = maxNMatrices;

// preallocate data structures
this->rows = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
this->cols = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );
this->offsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
this->rowOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
this->colOffsets = ( long * ) malloc( maxNMatrices * sizeof( long ) );
this->nnz = ( eslocal * ) malloc( maxNMatrices * sizeof( eslocal ) );

this->nMatrices = 0;
this->offsets[ 0 ] = 0;
this->rowOffsets[ 0 ] = 0;
this->colOffsets[ 0 ] = 0;
this->totalRows = 0;
this->totalCols = 0;
}

*/

/*SparseMatrixPack::SparseMatrixPack( const SparseMatrixPack& orig ) {

    this->device = orig.device;

    this->maxNMatrices = orig.maxNMatrices;
    this->nMatrices = orig.nMatrices;
    this->preallocSize = orig.preallocSize;

    // preallocate data structures
    this->matrix_values = orig.matrix_values;
    this->matrix_values_mic = orig.matrix_values_mic;
    this->rowInd = orig.rowInd;
    this->rowInd_mic = orig.rowInd_mic;
    this->colInd = orig.colInd;
    this->colInd_mic = orig.colInd_mic;
    this->rows = orig.rows;
    this->cols = orig.cols;
    this->nnz = orig.nnz;
    this->totalRows = orig.totalRows;
    this->totalCols = orig.totalCols;
    this->x_in_dim = orig.x_in_dim;
    this->y_out_dim = orig.y_out_dim;
    this->offsets = orig.offsets;
    this->rowOffsets = orig.rowOffsets;
    this->colOffsets = orig.colOffsets;
    this->nnz = orig.nnz;

    this->mic_x_in = orig.mic_x_in;
    this->mic_y_out = orig.mic_y_out;
    this->MICratio = orig.MICratio;
    this->elapsedTime = orig.elapsedTime;

    this->copiedToMIC = orig.copiedToMIC;
}*/

SparseMatrixPack::~SparseMatrixPack() {
    // free the MIC's memory
    if (this->copiedToMIC) {
#pragma offload_transfer target(mic:device) if(1) \ 
        nocopy( matrix_values_mic : alloc_if( 0 ) free_if( 1 ) targetptr ) \ 
       nocopy( rowInd_mic : alloc_if( 0 ) free_if( 1 ) targetptr ) \
        nocopy( colInd_mic : alloc_if( 0 ) free_if( 1 ) targetptr ) \
        nocopy( rows : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( cols : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( nnz : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( offsets : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( rowOffsets : alloc_if( 0 ) free_if(1) ) \
        nocopy( colOffsets : alloc_if( 0 ) free_if(1) ) \
       nocopy( mic_x_in : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( mic_y_out : alloc_if( 0 ) free_if( 1 ) ) \
        nocopy( elapsedTime : alloc_if( 0 ) free_if( 1 ) )
    }

    if (this->matrix_values != NULL )
        _mm_free( this->matrix_values );
    if (this->rows != NULL )
        delete [] this->rows ;
    if (this->cols != NULL )
        delete []  this->cols ;
    if (this->offsets != NULL )
        delete [] this->offsets;
    if (this->rowOffsets != NULL )
        delete [] this->rowOffsets;
    if (this->colOffsets != NULL )
        delete [] this->colOffsets;
    if (this->nnz != NULL )
        delete [] this->nnz;
    if ( this->mic_x_in != NULL )
        _mm_free( this->mic_x_in );
    if (this->mic_y_out != NULL )
        _mm_free( this->mic_y_out );
    if (this->elapsedTime != NULL) 
        delete [] elapsedTime;
    if ( this->rowInd != NULL ) 
        _mm_free( this->rowInd );
    if ( this->colInd != NULL )
        _mm_free( this->colInd );

}

void SparseMatrixPack::SetDevice(
        int device
        ) {
    this->device = device;
}

/*void SparseMatrixPack::PreparePack(
  eslocal i,
  eslocal nRows,
  eslocal nCols,
  eslocal nnz
  ) {

  this->rows[ i ] = nRows;
  this->cols[ i ] = nCols;
  this->rowOffsets[ i ] = this->totalRows;
  this->colOffsets[ i ] = this->totalCols;
  this->totalRows += nRows + 1;
  this->totalCols += nnz;

  this->nnz[ i ] = nnz;
  if ( i > 0 ) {
  this->offsets[ i ] = this->offsets[ i - 1 ] + this->nnz[ i - 1 ];
  }
  this->nMatrices++;
  }

  void SparseMatrixPack::AddSparseMatrix(
  eslocal i,
  double * matrixData
  ) {

  if ( i >= nMatrices ) {
  ESINFO(ERROR) << "Could not add matrix. Maximum number of matrices stored.";
  return;
  }

  eslocal size = this->nnz[i];




  if ( !this->packed[i] ) {
  memcpy( this->matrices + this->offsets[ i ],
  matrixData, sizeof( double ) * size );
  } else {
  for (long j = 0; j < this->cols[ i ]; j++ ) {
  long pos = ((double)( 1.0 +j ))*((double) j)/2.0;
  memcpy( this->matrices + this->offsets[ i ] + pos,
  matrixData + j * this->rows[ i ], sizeof( double ) * (j+1) );
  }
  }*/

void SparseMatrixPack::AddMatrices(
        SparseMatrix ** A, 
        eslocal n,
        eslocal device
        ) {

    this->device = device;

    this->rows = new eslocal[ n ];
    this->cols = new eslocal[ n ];
    this->offsets = new long[ n ];
    this->rowOffsets = new long[ n ];
    this->colOffsets = new long[ n ];
    this->x_in_offsets = new long[ n ];
    this->y_out_offsets = new long[ n ];
    this->nnz = new eslocal[ n ]; 

    this->offsets[ 0 ] = 0;
    this->x_in_offsets[ 0 ] = 0;
    this->y_out_offsets[ 0 ] = 0;
    this->totalRows = 0; 
    this->totalCols = 0;
    this->x_in_dim = 0;
    this->y_out_dim = 0;

    // iterate through matrices and sum up their sizes
    for ( eslocal i = 0; i < n ; ++i ) {
        this->rows[ i ] = A[i]->rows;
        this->cols[ i ] = A[i]->cols;
        this->rowOffsets[ i ] = this->totalRows;
        this->colOffsets[ i ] = this->totalCols;
        this->totalRows += (this->rows[i] + 1);
        this->totalCols += A[i]->nnz;
        this->x_in_dim += this->cols[i];
        this->y_out_dim += this->rows[i];

        this->nnz[ i ] = A[i]->nnz;
        this->preallocSize += A[i]->nnz;
        if ( i > 0 ) {
            this->offsets[ i ] = this->offsets[ i - 1 ] + this->nnz[ i - 1 ];
            this->x_in_offsets[ i ] = this->x_in_offsets[ i - 1 ] + this->cols[ i - 1 ];
            this->y_out_offsets[ i ] = this->y_out_offsets[ i - 1 ] + this->rows[ i - 1 ];
        }
        this->nMatrices++;
    }

    this->matrix_values = ( double * ) _mm_malloc( this->preallocSize * sizeof( double ), 64 ); 
    this->rowInd = ( eslocal * ) _mm_malloc( this->totalRows * sizeof( eslocal ), 64 );
    this->colInd = ( eslocal * ) _mm_malloc( this->totalCols * sizeof( eslocal ), 64 );


    for ( eslocal i = 0 ; i < n ; ++i ) {
        memcpy( this->matrix_values + this->offsets[ i ], &(A[i]->CSR_V_values[0]), A[i]->nnz * sizeof( double ) );
        memcpy( this->rowInd + this->rowOffsets[ i ], &(A[i]->CSR_I_row_indices[0]), (A[i]->rows + 1) * sizeof( eslocal ) );   
        memcpy( this->colInd + this->colOffsets[ i ], &(A[i]->CSR_J_col_indices[0]), A[i]->nnz * sizeof( eslocal ) );   
    }
}


void SparseMatrixPack::SetX(
        long vector,
        long position,
        double value
        ) {
    if ( this->mic_x_in != NULL ) {
        this->mic_x_in[position + this->x_in_offsets[ vector ]] = value;
    } else {
        ESINFO(ERROR) << "Could not set vector element. Vector not yet allocated.";
    }
}

void SparseMatrixPack::GetY(
        long vector,
        std::SEQ_VECTOR <double> & y
        ) {
    if ( this->mic_y_out != NULL ) {
        memcpy( &(y[0]), this->mic_y_out + this->y_out_offsets[ vector ],
                this->rows[vector] * sizeof( double ) );
    } else {
        ESINFO(ERROR) << "Could not copy vector. Vector not yet allocated.";
    }
}


void SparseMatrixPack::CopyToMIC( ) {
    double * tmp_val;
    eslocal * tmp_rowInd;
    eslocal * tmp_colInd;
    // allocate targetptr array on MIC
#pragma offload_transfer target(mic:device) \
    nocopy(tmp_val : length( this->preallocSize ) alloc_if(1) free_if(0) targetptr) \
    nocopy(tmp_rowInd : length( this->totalRows ) alloc_if(1) free_if(0) targetptr) \
    nocopy(tmp_colInd : length( this->totalCols ) alloc_if(1) free_if(0) targetptr) \
    in(this : alloc_if(1) free_if(0))

    // preallocated input/output buffers
    //mic_x_in = ( double * ) malloc( totalCols * sizeof( double ) );
    //mic_y_out = ( double * ) malloc( totalRows * sizeof( double ) );
    this->matrix_values_mic = tmp_val;
    this->rowInd_mic = tmp_rowInd;
    this->colInd_mic = tmp_colInd;

#pragma offload_transfer target(mic:device) if(1) \
    in( matrix_values : length( this->preallocSize ) into( this->matrix_values_mic ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rowInd : length( this->totalRows ) into( this->rowInd_mic ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( colInd : length( this->totalCols ) into( this->colInd_mic ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rows : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
    in( cols : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
    in( nnz : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
    in( offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
    in( rowOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
    in( colOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
    in( x_in_offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
    in( y_out_offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 )) \
    in( mic_x_in : length( this->x_in_dim ) alloc_if( 1 ) free_if( 0 ) ) \
    in( mic_y_out : length( this->y_out_dim ) alloc_if( 1 ) free_if( 0 ) ) \
    in( elapsedTime : length( 1 ) alloc_if( 1 ) free_if( 0 ) ) \
    in( this : length(0) alloc_if( 0 ) free_if( 0 ) )

    this->copiedToMIC = true;
    if ( !configuration.load_balancing ) {
    //    _mm_free( this->matrix_values );
    //    _mm_free( this->rowInd );
    //    _mm_free( this->colInd );
    //    this->matrix_values = NULL;
    //    this->rowInd = NULL; 
    //    this->colInd = NULL;
    }
}

void SparseMatrixPack::FactorizeMIC( ) {


    pt = new void**[nMatrices];
    iparm = new MKL_INT*[nMatrices];
    dparm = new double*[nMatrices];
    perm = new MKL_INT*[nMatrices];
    error = new MKL_INT[nMatrices];

    for (eslocal i = 0; i < nMatrices; i++) {

        pt[i] = new void*[64];
        iparm[i] = new MKL_INT[64];
        dparm[i] = new double[65];
        perm[i] = new MKL_INT[ rows[i] ];
        for (eslocal j = 0; j < 64; j++) {
            iparm[i][j]=0;
        }
        for (eslocal j = 0; j < rows[i] ; j++) {
            perm[i][j] = j+1; 
        }
        for (eslocal j = 0 ;j < 64; j++) {
            pt[i][j] = 0;
        }
    }

#pragma offload_transfer target(mic:device) \
    in(pt : length(nMatrices) alloc_if(1) free_if(0)) \
    in(iparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(dparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(perm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(error : length(nMatrices) alloc_if(1) free_if(0)) \
    in(this : length(0) alloc_if(0) free_if(0))
    for (eslocal i = 0; i<nMatrices; i++) {
        void **ptPointer = pt[i];
        MKL_INT* iparmPointer = iparm[i];
        double* dparmPointer = dparm[i];
        MKL_INT *permPointer = perm[i];
#pragma offload target(mic:device) \
        in( ptPointer : length(64) alloc_if(1) free_if(0)) \
        in( iparmPointer : length(64) alloc_if(1) free_if(0)) \
        in( dparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( permPointer : length(rows[i]) alloc_if(1) free_if(0)) \
        in(iparm : length(0) alloc_if(0) free_if(0)) \
        in(dparm : length(0) alloc_if(0) free_if(0)) \
        in(perm : length(0) alloc_if(0) free_if(0)) \
        in(pt : length(0) alloc_if(0) free_if(0)) \
        in(this : length(0) alloc_if(0) free_if(0))
        {
            pt[i] = ptPointer;
            iparm[i] = iparmPointer;
            dparm[i] = dparmPointer;
            perm[i] = permPointer;
        }
    }


#pragma offload target( mic : device ) \
    in( matrix_values_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rowInd_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( colInd_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( nnz : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( colOffsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( pt : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( iparm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( dparm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( error : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( perm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) 
    {
        MKL_INT phase = 11;
        MKL_INT idum;
        MKL_INT mnum = 1;
        MKL_INT maxfct = 1;
        double ddum; 
        MKL_INT m_nRhs = 0;
        MKL_INT msglvl =  0;
        MKL_INT mtype = 2;


#pragma omp parallel
        {
#pragma omp for 
            for (eslocal i = 0 ; i < nMatrices ; ++i ) {
                iparm[i][2] = 0;

                iparm[i][1-1] = 1;         /* No solver default */
                iparm[i][2-1] = 2;         /* Fill-in reordering from METIS */
                iparm[i][5-1] = 0;         /* No user fill-in reducing permutation */
                iparm[i][6-1] = 0;         /* Write solution to x */
                iparm[i][10-1] = 8;   /* Perturb the pivot elements with 1E-13 */
                iparm[i][11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
                iparm[i][13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
                iparm[i][14-1] = 0;        /* Output: Number of perturbed pivots */
                iparm[i][18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
                iparm[i][19-1] = -1;       /* Output: Mflops for LU factorization */
                iparm[i][24-1] = 0;        /* Parallel factorization control */
                iparm[i][25-1] = 0;        /* Parallel forward/backward solve control */
                iparm[i][27-1] = 0;
                iparm[i][28-1] = 0;
                iparm[i][36-1] = 0;        /* Use Schur complement */


                pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
                        &rows[i], matrix_values_mic + offsets[i], 
                        rowInd_mic + rowOffsets[i], colInd_mic + colOffsets[i], perm[i], 
                        &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
            }
            bool initialized = true;
            for (eslocal i=0; i < nMatrices; i++) {
                if (error[i] != 0) {
                    initialized = false;
                    std::cerr << "ERROR during symbolic factorization of matrix " << i << " : " << error[i] << "\n";
                }
            }
            if (!initialized) {
                exit(EXIT_FAILURE);
            }

#ifdef DEBUG
            preslocalf ("\nReordering completed ... ");
#endif

            /* -------------------------------------------------------------------- */
            /* .. Numerical factorization. */
            /* -------------------------------------------------------------------- */
#pragma omp single
            {
                phase = 22;
            }
#pragma omp for
            for (eslocal i = 0; i < nMatrices; i++) {
                pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
                        &rows[i], matrix_values_mic + offsets[i], 
                        rowInd_mic + rowOffsets[i], colInd_mic + colOffsets[i], perm[i], 
                        &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
            }
            bool m_factorized = true;
            for (eslocal i=0; i < nMatrices; i++) {
                if (error[i] != 0) {
                    m_factorized = false;
                    std::cerr << "ERROR during numeric factorization of matrix " << i << "\n";
                }
            }
            if (m_factorized!=true) {
                exit(EXIT_FAILURE);
            }

#ifdef DEBUG
            preslocalf ("\nFactorization completed ... ");
#endif

        }
    }
}

void SparseMatrixPack::SolveMIC_Start( 
        ) {
#pragma offload target( mic : device ) \
    in( matrix_values_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rowInd_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( colInd_mic : length( 0 ) alloc_if( 0 ) free_if( 0 ) targetptr ) \
    in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( nnz : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( colOffsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( x_in_offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( y_out_offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( pt : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( iparm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( dparm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( error : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( perm : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( mic_x_in : length( x_in_dim ) alloc_if( 0 ) free_if( 0 ) ) \
    in( mic_y_out : length( 0 ) alloc_if(0 ) free_if( 0 ) ) \
    in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    signal( mic_y_out )
    {

        double ddum     = 0;
        MKL_INT idum    = 0; 
        MKL_INT n_rhs   = 1;
        MKL_INT msglvl  = 0;
        MKL_INT mnum    = 1;
        MKL_INT maxfct  = 1;
        MKL_INT mtype   = 2;

        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */

#pragma omp parallel 
        {
            eslocal myPhase = 33;
#pragma omp for 
            for (eslocal i = 0; i < nMatrices; i++) {
                MKL_INT ip5backup = iparm[i][5];
                iparm[i][ 6 - 1 ] = 1;
               double * out =  mic_y_out + y_out_offsets[i]; 
                pardiso ( pt[i], &maxfct, &mnum, &mtype, &myPhase,
                        &rows[i], matrix_values_mic + offsets[i],
                        rowInd_mic + rowOffsets[i], colInd_mic + colOffsets[i], 
                        &idum, &n_rhs, iparm[i], &msglvl, mic_x_in + x_in_offsets[i], 
                        out, &error[i] );
                iparm[i][5] = ip5backup;
            }
        }
        bool err = false;
        for (eslocal i = 0; i < nMatrices; i++) {
            if (error[i]!=0) {
                err = true;
            }
        }
        if (err)
        {
            exit (3);
        }
    }
}

void SparseMatrixPack::SolveMIC_Sync( 
        ) {
#pragma offload_wait target( mic : device ) wait( mic_y_out ) 
#pragma offload_transfer target( mic : device ) \
    out( mic_x_in : into(mic_y_out) length( y_out_dim ) alloc_if( 0 ) free_if( 0 ) )
}


}
