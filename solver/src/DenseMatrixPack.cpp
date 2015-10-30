#include "DenseMatrixPack.h"

DenseMatrixPack::DenseMatrixPack() {

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

}

DenseMatrixPack::DenseMatrixPack(
  long maxNMatrices,
  long preallocSize,
  int device
) {

  this->device = device;

  this->maxNMatrices = maxNMatrices;
  this->preallocSize = preallocSize;

  // preallocate data structures
  this->matrices = ( double * ) malloc( preallocSize * sizeof( double ) );
  this->rows = ( int * ) malloc( maxNMatrices * sizeof( int ) );
  this->cols = ( int * ) malloc( maxNMatrices * sizeof( int ) );
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
  this->rows = ( int * ) malloc( maxNMatrices * sizeof( int ) );
  this->cols = ( int * ) malloc( maxNMatrices * sizeof( int ) );
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

  this->device = orig.device;

  this->maxNMatrices = orig.maxNMatrices;
  this->preallocSize = orig.preallocSize;

  // preallocate data structures
  this->matrices = orig.matrices;
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
}

DenseMatrixPack::~DenseMatrixPack() {
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

void DenseMatrixPack::PreparePack(
  int i,
  int nRows,
  int nCols,
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
    printf( "Could not add matrix. Maximum number of matrices stored." );
  } else if ( size > freeSpace ) {
    printf( "Could not add matrix. Not enough allocated memory." );
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
  int i,
  double * matrixData
) {

  if ( i >= nMatrices ) {
    std::cout << "Could not add matrix. Maximum number of matrices stored."
      << std::endl;
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
    std::cout << "Could not set vector element. Vector not yet allocated."
      << std::endl;
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
    std::cout << "Could not copy vector. Vector not yet allocated."
      << std::endl;
  }
}

void DenseMatrixPack::CopyToMIC( ) {
#ifdef MIC
  long matrixLength = preallocSize - freeSpace;

  // preallocated input/output buffers
  //mic_x_in = ( double * ) malloc( totalCols * sizeof( double ) );
  //mic_y_out = ( double * ) malloc( totalRows * sizeof( double ) );
#pragma offload_transfer target(mic:device) if(1) \
  in( matrices : length( matrixLength ) alloc_if( 1 ) free_if( 0 ) ) \
  in( rows : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
  in( cols : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
  in( offsets : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
  in( rowOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
  in( colOffsets : length( nMatrices ) alloc_if( 1 ) free_if(0) ) \
  in( lengths : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
  in( packed : length( nMatrices ) alloc_if( 1 ) free_if( 0 ) ) \
  in( mic_x_in : length( totalCols ) alloc_if( 1 ) free_if( 0 ) ) \
  in( mic_y_out : length( totalRows ) alloc_if( 1 ) free_if( 0 ) ) \
  in( this : alloc_if( 1 ) free_if( 0 ) )

#endif
}

void DenseMatrixPack::DenseMatsVecs(
  char T_for_transpose_N_for_not_transpose
) {

			double alpha = 1.0;
			double beta  = 0.0;
			int one = 1;

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

void DenseMatrixPack::DenseMatsVecsMIC(
  char T_for_transpose_N_for_not_transpose
) {
#ifdef MIC
	#pragma offload target(mic:device) if(1) \
		in( matrices : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
		in( mic_x_in : length( totalCols ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
		out( mic_y_out : length( totalRows ) alloc_if( 0 ) free_if( 0 ) ) \
    in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
		{
			double alpha = 1.0;
			double beta  = 0.0;
			int one = 1;

    //  omp_set_num_threads(24);

      #pragma omp parallel for //schedule(dynamic,10)
      for ( int i = 0 ; i < nMatrices; i++ ) {

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

#endif

}

void DenseMatrixPack::DenseMatsVecsMIC_Start(
  char T_for_transpose_N_for_not_transpose
) {
#ifdef MIC
	#pragma offload target(mic:device) if(1) \
		in( matrices : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
		in( mic_x_in : length( totalCols ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( mic_y_out : length( 0 ) alloc_if( 0 ) free_if( 0 )  ) \
		signal( mic_y_out )\
    in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
		{
			double alpha = 1.0;
			double beta  = 0.0;
			int one = 1;

    //  omp_set_num_threads(24);

      #pragma omp parallel for schedule(dynamic)
      for ( int i = 0 ; i < nMatrices; i++ ) {

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
#endif
}

void DenseMatrixPack::DenseMatsVecsMIC_Sync( ) {
#ifdef MIC
  #pragma offload target(mic:device) if(1) \
    in( matrices : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( mic_x_in : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rows : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( cols : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( offsets : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( rowOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( colOffsets : length( 0 ) alloc_if( 0 ) free_if(0) ) \
    in( lengths : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    in( packed : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
    out( mic_y_out : length( totalCols ) alloc_if( 0 ) free_if( 0 )  ) \
    wait( mic_y_out )\
    in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
    {}
#endif
}
