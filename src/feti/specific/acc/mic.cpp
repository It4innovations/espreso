#include "mic.h"

using namespace espreso;

SparseSolverMIC::SparseSolverMIC() {
    device = 0;
    isOffloaded = false;
    keep_factors = true;
    initialized = false;
    USE_FLOAT = false;

    nMatrices = 0;

    mtype = 2;
    tmp_sol_fl1 = 0;
    tmp_sol_fl1 = 0;

}

SparseSolverMIC::~SparseSolverMIC() {
    this->Clear();    
}

void SparseSolverMIC::Clear() {
    if (isOffloaded) {
        for (esint i = 0; i<nMatrices; i++) {
            double* valPointer = CSR_V_values[i];
            MKL_INT* iPointer = CSR_I_row_indices[i];
            MKL_INT* jPointer = CSR_J_col_indices[i];
            void **ptPointer = pt[i];
            MKL_INT* iparmPointer = iparm[i];
            double* dparmPointer = dparm[i];
            MKL_INT *permPointer = perm[i];
#pragma offload_transfer target(mic:device) \
            in( valPointer : length(0) alloc_if(0) free_if(1)) \
            in( iPointer : length(0) alloc_if(0) free_if(1)) \
            in( jPointer : length(0) alloc_if(0) free_if(1)) \
            in( ptPointer : length(0) alloc_if(0) free_if(1)) \
            in( iparmPointer : length(0) alloc_if(0) free_if(1)) \
            in( dparmPointer : length(0) alloc_if(0) free_if(1)) \
            in( permPointer : length(0) alloc_if(0) free_if(1))
        }


#pragma offload_transfer target(mic:device) \
        in(rows : length(0) alloc_if(0) free_if(1) ) \
        in(cols : length(0) alloc_if(0) free_if(1) ) \
        in(nnz : length(0) alloc_if(0) free_if(1) ) \
        in(CSR_I_row_indices_size : length(0) alloc_if(0) free_if(1)) \
        in(CSR_J_col_indices_size : length(0) alloc_if(0) free_if(1)) \
        in(CSR_V_values_size : length(0) alloc_if(0) free_if(1)) \
        in(pt : length(0) alloc_if(0) free_if(1)) \
        in(iparm : length(0) alloc_if(0) free_if(1)) \
        in(dparm : length(0) alloc_if(0) free_if(1)) \
        in(perm : length(0) alloc_if(0) free_if(1)) \
        in(error : length(0) alloc_if(0) free_if(1)) \
        in(m_Kplus_size : length(0) alloc_if(0) free_if(1))

        _mm_free(rows);
        _mm_free(cols);
        _mm_free(nnz);
        _mm_free(CSR_I_row_indices_size);
        _mm_free(CSR_J_col_indices_size);
        _mm_free(CSR_V_values_size);

        delete [] m_Kplus_size;
        delete [] error;
        for (esint i = 0 ; i < nMatrices; i++) {
            if(import_with_copy) {
                _mm_free(CSR_I_row_indices[i]);
                _mm_free(CSR_J_col_indices[i]);
                _mm_free(CSR_V_values[i]);
            }
            delete [] pt[i];
            delete [] iparm[i];
            delete [] dparm[i];
            _mm_free(perm[i]);
        }
        delete [] pt;
        delete [] iparm;
        delete [] dparm;
        delete [] perm;
    }


}


void SparseSolverMIC::ImportMatrices_wo_Copy(SparseMatrix ** A, esint nMatrices, esint mic) {

    USE_FLOAT = false;
    this->nMatrices = nMatrices;

    this->device = mic;
    this->A = A;

/*
    rows = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    cols = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    nnz =  (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_J_col_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_V_values_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices = new MKL_INT*[nMatrices];
    CSR_J_col_indices = new MKL_INT*[nMatrices];
    CSR_V_values = new double*[nMatrices];

    pt = new void**[nMatrices];
    iparm = new MKL_INT*[nMatrices];
    dparm = new double*[nMatrices];
    perm = new MKL_INT*[nMatrices];
    error = new MKL_INT[nMatrices];
    m_factorized = 0;
    m_Kplus_size = new MKL_INT[nMatrices];

    for (esint i = 0; i < nMatrices; i++) {
        rows[i] = A[i]->rows;
        cols[i] = A[i]->cols;
        nnz[i] = A[i]->nnz;
        m_Kplus_size[i] = A[i]->rows;
        CSR_I_row_indices_size[i] = A[i]->CSR_I_row_indices.size();
        CSR_J_col_indices_size[i] = A[i]->CSR_J_col_indices.size();
        CSR_V_values_size[i] = A[i]->CSR_V_values.size();

        CSR_I_row_indices[i] = &(A[i]->CSR_I_row_indices[0]);
        CSR_J_col_indices[i] = &(A[i]->CSR_J_col_indices[0]);
        CSR_V_values[i] = &(A[i]->CSR_V_values[0]);

        pt[i] = new void*[64];
        iparm[i] = new MKL_INT[65];
        dparm[i] = new double[65];
        perm[i] = (MKL_INT*) _mm_malloc(rows[i] * sizeof(MKL_INT), 64);
        for (esint j = 0; j < 65; j++) {
            iparm[i][j]=0;
        }
        for (esint j = 0; j < rows[i] ; j++) {
            perm[i][j] = 0; 
        }
        for (esint j = 0 ;j < 64; j++) {
            pt[i][j] = 0;
        }
    }
    // send data to MIC

#pragma offload_transfer target(mic:device) \
    in(rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(cols : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(nnz : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(CSR_I_row_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_J_col_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_V_values_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_V_values : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_I_row_indices : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_J_col_indices : length(nMatrices) alloc_if(1) free_if(0)) \
    in(pt : length(nMatrices) alloc_if(1) free_if(0)) \
    in(iparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(dparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(perm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(error : length(nMatrices) alloc_if(1) free_if(0)) \
    in(m_Kplus_size : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(this : alloc_if(1) free_if(0))
    for (esint i = 0; i<nMatrices; i++) {
        double* valPointer = CSR_V_values[i];
        MKL_INT* iPointer = CSR_I_row_indices[i];
        MKL_INT* jPointer = CSR_J_col_indices[i];
        void **ptPointer = pt[i];
        MKL_INT* iparmPointer = iparm[i];
        double* dparmPointer = dparm[i];
        MKL_INT *permPointer = perm[i];
#pragma offload target(mic:device) \
        in( valPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( iPointer : length(rows[i]+1) alloc_if(1) free_if(0)) \
        in( jPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( ptPointer : length(64) alloc_if(1) free_if(0)) \
        in( iparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( dparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( permPointer : length(rows[i]) alloc_if(1) free_if(0)) \
        in(CSR_V_values : length(0) alloc_if(0) free_if(0)) \
        in(CSR_I_row_indices : length(0) alloc_if(0) free_if(0)) \
        in(CSR_J_col_indices : length(0) alloc_if(0) free_if(0)) \
        in(iparm : length(0) alloc_if(0) free_if(0)) \
        in(dparm : length(0) alloc_if(0) free_if(0)) \
        in(perm : length(0) alloc_if(0) free_if(0)) \
        in(pt : length(0) alloc_if(0) free_if(0)) \
        in(this : length(0) alloc_if(0) free_if(0))
        {
            CSR_V_values[i] = valPointer;
            CSR_I_row_indices[i] = iPointer;
            CSR_J_col_indices[i] = jPointer;
            pt[i] = ptPointer;
            iparm[i] = iparmPointer;
            dparm[i] = dparmPointer;
            perm[i] = permPointer;
        }
    }
    isOffloaded = true;
*/
}




void SparseSolverMIC::ImportMatrices(SparseMatrix ** A, esint nMatrices, esint mic) {

    USE_FLOAT = false;
    this->nMatrices = nMatrices;

    this->device = mic;
    this->A = A;
/*
    rows = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    cols = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    nnz =  (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_J_col_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_V_values_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices = new MKL_INT*[nMatrices];
    CSR_J_col_indices = new MKL_INT*[nMatrices];
    CSR_V_values = new double*[nMatrices];

    pt = new void**[nMatrices];
    iparm = new MKL_INT*[nMatrices];
    dparm = new double*[nMatrices];
    perm = new MKL_INT*[nMatrices];
    error = new MKL_INT[nMatrices];
    m_factorized = 0;
    m_Kplus_size = new MKL_INT[nMatrices];

    for (esint i = 0; i < nMatrices; i++) {
        rows[i] = A[i]->rows;
        cols[i] = A[i]->cols;
        nnz[i] = A[i]->nnz;
        m_Kplus_size[i] = A[i]->rows;
        CSR_I_row_indices_size[i] = A[i]->CSR_I_row_indices.size();
        CSR_J_col_indices_size[i] = A[i]->CSR_J_col_indices.size();
        CSR_V_values_size[i] = A[i]->CSR_V_values.size();

        CSR_I_row_indices[i] = (MKL_INT*) _mm_malloc(CSR_I_row_indices_size[i] * sizeof(MKL_INT), 64);
        CSR_J_col_indices[i] = (MKL_INT*) _mm_malloc(CSR_J_col_indices_size[i] * sizeof(MKL_INT), 64);
        CSR_V_values[i] = (double*) _mm_malloc(CSR_V_values_size[i] * sizeof(double), 64);

        copy(A[i]->CSR_I_row_indices.begin(), A[i]->CSR_I_row_indices.end(), CSR_I_row_indices[i]);
        copy(A[i]->CSR_J_col_indices.begin(), A[i]->CSR_J_col_indices.end(), CSR_J_col_indices[i]);
        copy(A[i]->CSR_V_values.begin(), A[i]->CSR_V_values.end(), CSR_V_values[i]);

        pt[i] = new void*[65];
        iparm[i] = new MKL_INT[65];
        dparm[i] = new double[65];
        perm[i] = (MKL_INT*) _mm_malloc(rows[i] * sizeof(MKL_INT), 64);
        for (esint j = 0; j < 65; j++) iparm[j]=0;
    }
    import_with_copy = true;

    // send data to MIC

#pragma offload_transfer target(mic:device) \
    in(rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(cols : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(nnz : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(CSR_I_row_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_J_col_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_V_values_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(pt : length(nMatrices) alloc_if(1) free_if(0)) \
    in(iparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(dparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(perm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(error : length(nMatrices) alloc_if(1) free_if(0)) \
    in(m_Kplus_size : length(nMatrices) alloc_if(1) free_if(0))

    for (esint i = 0; i<nMatrices; i++) {
        double* valPointer = CSR_V_values[i];
        MKL_INT* iPointer = CSR_I_row_indices[i];
        MKL_INT* jPointer = CSR_J_col_indices[i];
        void **ptPointer = pt[i];
        MKL_INT* iparmPointer = iparm[i];
        double* dparmPointer = dparm[i];
        MKL_INT *permPointer = perm[i];
#pragma offload target(mic:device) \
        in( valPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( iPointer : length(rows[i]+1) alloc_if(1) free_if(0)) \
        in( jPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( ptPointer : length(64) alloc_if(1) free_if(0)) \
        in( iparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( dparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( permPointer : length(rows[i]) alloc_if(1) free_if(0))
        {
            CSR_V_values[i] = valPointer;
            CSR_I_row_indices[i] = iPointer;
            CSR_J_col_indices[i] = jPointer;
            pt[i] = ptPointer;
            iparm[i] = iparmPointer;
            dparm[i] = dparmPointer;
            perm[i] = permPointer;
        }
    }
    isOffloaded = true;
*/
}

void SparseSolverMIC::ImportMatrices_fl(SparseMatrix ** A, esint nMatrices, esint mic) {

    USE_FLOAT = true;

    this->device = mic;
    this->nMatrices = nMatrices;
    this->A = A;

/*    rows = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    cols = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    nnz =  (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_J_col_indices_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_V_values_size = (MKL_INT*) _mm_malloc(nMatrices * sizeof(MKL_INT), 64);
    CSR_I_row_indices = new MKL_INT*[nMatrices];
    CSR_J_col_indices = new MKL_INT*[nMatrices];
    CSR_V_values_fl = new float*[nMatrices];

    pt = new void**[nMatrices];
    iparm = new MKL_INT*[nMatrices];
    dparm = new double*[nMatrices];
    error = new MKL_INT[nMatrices];
    m_factorized = 0;
    m_Kplus_size = new MKL_INT[nMatrices];


    for (esint i = 0; i < nMatrices; i++) {
        rows[i] = A[i]->rows;
        cols[i] = A[i]->cols;
        nnz[i] = A[i]->nnz;
        m_Kplus_size[i] = A[i]->rows;
        CSR_I_row_indices_size[i] = A[i]->CSR_I_row_indices.size();
        CSR_J_col_indices_size[i] = A[i]->CSR_J_col_indices.size();
        CSR_V_values_size[i] = A[i]->CSR_V_values.size();

        CSR_I_row_indices[i] = (MKL_INT*) _mm_malloc(CSR_I_row_indices_size[i] * sizeof(MKL_INT), 64);
        CSR_J_col_indices[i] = (MKL_INT*) _mm_malloc(CSR_J_col_indices_size[i] * sizeof(MKL_INT), 64);
        CSR_V_values_fl[i] = (float*) _mm_malloc(CSR_V_values_size[i] * sizeof(float), 64);

        copy(A[i]->CSR_I_row_indices.begin(), A[i]->CSR_I_row_indices.end(), CSR_I_row_indices[i]);
        copy(A[i]->CSR_J_col_indices.begin(), A[i]->CSR_J_col_indices.end(), CSR_J_col_indices[i]);
        copy(A[i]->CSR_V_values.begin(), A[i]->CSR_V_values.end(), CSR_V_values_fl[i]);

        pt[i] = new void*[65];
        iparm[i] = new MKL_INT[65];
        dparm[i] = new double[65];

    }
    import_with_copy = true;

    //send data to MIC

#pragma offload_transfer target(mic:device) \
    in(rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(cols : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(nnz : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(CSR_I_row_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_J_col_indices_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(CSR_V_values_size : length(nMatrices) alloc_if(1) free_if(0)) \
    in(pt : length(nMatrices) alloc_if(1) free_if(0)) \
    in(iparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(dparm : length(nMatrices) alloc_if(1) free_if(0)) \
    in(error : length(nMatrices) alloc_if(1) free_if(0)) \
    in(m_Kplus_size : length(nMatrices) alloc_if(1) free_if(0))


    for (esint i = 0; i<nMatrices; i++) {
        double* valPointer = CSR_V_values[i];
        MKL_INT* iPointer = CSR_I_row_indices[i];
        MKL_INT* jPointer = CSR_J_col_indices[i];
        void **ptPointer = pt[i];
        MKL_INT* iparmPointer = iparm[i];
        double* dparmPointer = dparm[i];


#pragma offload target(mic:device) \
        in( valPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( iPointer : length(rows[i]+1) alloc_if(1) free_if(0)) \
        in( jPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( ptPointer : length(64) alloc_if(1) free_if(0)) \
        in( iparmPointer : length(65) alloc_if(1) free_if(0)) \
        in( dparmPointer : length(65) alloc_if(1) free_if(0))


        {
            CSR_V_values[i] = valPointer;
            CSR_I_row_indices[i] = iPointer;
            CSR_J_col_indices[i] = jPointer;
            pt[i] = ptPointer;
            iparm[i] = iparmPointer;
            dparm[i] = dparmPointer;

        }
    }

    isOffloaded = true;
    */
}

void SparseSolverMIC::Factorization(const std::string &str) {
    // reordering and symbolic factorization of all matrices at MIC in parallel
    //

#pragma offload target(mic:device) \
    in(rows : length(0) alloc_if(0) free_if(0) ) \
    in(cols : length(0) alloc_if(0) free_if(0) ) \
    in(nnz : length(0) alloc_if(0) free_if(0) ) \
    in(CSR_I_row_indices : length(0) alloc_if(0) free_if(0)) \
    in(CSR_J_col_indices : length(0) alloc_if(0) free_if(0)) \
    in(CSR_V_values : length(0) alloc_if(0) free_if(0)) \
    in(pt : length(0) alloc_if(0) free_if(0)) \
    in(iparm : length(0) alloc_if(0) free_if(0)) \
    in(dparm : length(0) alloc_if(0) free_if(0)) \
    in(error : length(0) alloc_if(0) free_if(0)) \ 
    in(m_Kplus_size : length(0) alloc_if(0) free_if(0)) \
        in(perm : length(0) alloc_if(0) free_if(0)) \
        in(this : length(0) alloc_if(0) free_if(0))
        {
            phase = 11;
            MKL_INT idum;
            mnum = 1;
            maxfct = 1;
            double ddum;
            m_nRhs = 0;
            msglvl = 0;
#pragma omp parallel
            {
#pragma omp for
                for (esint i = 0; i < nMatrices; i++) {
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
                    iparm[i][28-1] = 0;
                    iparm[i][36-1] = 0;        /* Use Schur complement */

                    if (USE_FLOAT) {
                        iparm[i][27] = 1; // run PARDISO in float
                        pardiso(pt[i], &maxfct, &mnum, &mtype, &phase,
                                &rows[i], CSR_V_values_fl[i], CSR_I_row_indices[i], CSR_J_col_indices[i], perm[i], &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
                    } else {
                        pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
                                &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i], perm[i], &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
                    }

                }
                initialized = true;            
                for (esint i=0; i < nMatrices; i++) {
                    if (error[i] != 0) {
                        initialized = false;
                        std::cerr << "ERROR during symbolic factorization of matrix " << i <<": " << str << "\n";
                    }
                }
                if (!initialized) {
                    exit(EXIT_FAILURE);
                }

#ifdef DEBUG
                presintf ("\nReordering completed ... ");
#endif

                /* -------------------------------------------------------------------- */
                /* .. Numerical factorization. */
                /* -------------------------------------------------------------------- */
#pragma omp single
                {
                    phase = 22;
                }
#pragma omp for
                for (esint i = 0; i < nMatrices; i++) {
                    if (USE_FLOAT) {
                        pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
                                &rows[i], CSR_V_values_fl[i], CSR_I_row_indices[i], CSR_J_col_indices[i], perm[i], &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
                    } else {
                        pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
                                &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i], perm[i], &m_nRhs, iparm[i], &msglvl, &ddum, &ddum, &error[i]);
                    }
                    m_factorized = 1;
                    for (esint i=0; i < nMatrices; i++) {
                        if (error[i] != 0) {
                            m_factorized = 0;
                            std::cerr << "ERROR during numeric factorization of matrix " << i <<": " << str << "\n";
                        }
                    }
                    if (m_factorized!=1) {
                        exit(EXIT_FAILURE);
                    }

#ifdef DEBUG
                    presintf ("\nFactorization completed ... ");
#endif
                }
            }
        }
    if (USE_FLOAT) {
        tmp_sol_fl1 = new float*[nMatrices];
        tmp_sol_fl2 = new float*[nMatrices];
        for (esint i=0; i < nMatrices; i++) {
            tmp_sol_fl1[i] = (float*) _mm_malloc(m_Kplus_size[i]*sizeof(float), 64);
            tmp_sol_fl2[i] = (float*) _mm_malloc(m_Kplus_size[i]*sizeof(float), 64);
        }
        // send data persistently to mic so it does not have to be allocated
        // everytime the solve is needed
#pragma offload_transfer target(mic:device) \
        in(tmp_sol_fl1 : length(nMatrices) alloc_if(1) free_if(0)) \
        in(tmp_sol_fl2 : length(nMatrices) alloc_if(1) free_if(0))
        for (esint i = 0; i<nMatrices; i++) {
            float * tmp1 = tmp_sol_fl1[i];
            float * tmp2 = tmp_sol_fl2[i];
            esint tmpLength = m_Kplus_size[i];
#pragma offload target(mic:device) \
            in(tmp1 : length(tmpLength) alloc_if(1) free_if(0)) \
            in(tmp2 : length(tmpLength) alloc_if(1) free_if(0)) \
            in(this : length(0) alloc_if(0) free_if(0))
            {
                tmp_sol_fl1[i] = tmp1;
                tmp_sol_fl2[i] = tmp2;
            }
        }
    } else {
        tmp_sol_d1 = new double*[nMatrices];
        tmp_sol_d2 = new double*[nMatrices];
        for (esint i=0; i<nMatrices; i++) {
            tmp_sol_d1[i] = (double*) _mm_malloc(m_Kplus_size[i]*sizeof(double),64);
            tmp_sol_d2[i] = (double*) _mm_malloc(m_Kplus_size[i]*sizeof(double),64);
        }
        // send data persistently to mic so it does not have to be allocated
        // everytime the solve is needed
#pragma offload_transfer target(mic:device) \
        in(tmp_sol_d1 : length(nMatrices) alloc_if(1) free_if(0)) \
        in(tmp_sol_d2 : length(nMatrices) alloc_if(1) free_if(0)) \
        in(this : length(0) alloc_if(0) free_if(0))
        for (esint i = 0; i<nMatrices; i++) {
            double * tmp1 = tmp_sol_d1[i];
            double * tmp2 = tmp_sol_d2[i];
            esint tmpLength = m_Kplus_size[i];
#pragma offload target(mic:device) \
            in(tmp_sol_d1 : length(0) alloc_if(0) free_if(0)) \
            in(tmp_sol_d2 : length(0) alloc_if(0) free_if(0)) \
            in(tmp1 : length(tmpLength) alloc_if(1) free_if(0)) \
            in(tmp2 : length(tmpLength) alloc_if(1) free_if(0)) \
            in(this : length(0) alloc_if(0) free_if(0))
            {
                tmp_sol_d1[i] = tmp1;
                tmp_sol_d2[i] = tmp2;
            }
        }

    }
    
    int totalLength = 0;
    for (esint i = 0 ; i < nMatrices; i++) {
        totalLength += m_Kplus_size[i];
    }
    vectors = (double*) _mm_malloc(totalLength * sizeof(double), 64);
    vectors_out = (double*) _mm_malloc(totalLength * sizeof(double), 64);
#pragma offload_transfer target(mic:device) \
    in(vectors : length(totalLength) alloc_if(1) free_if(0)) \ 
    in(vectors_out : length(totalLength) alloc_if(1) free_if(0)) 

    m_factorized = true;
}


void SparseSolverMIC::Factorization(
    const std::string &str,
    SparseMatrixPack &factors_out
    ) {

    if ( this->A != NULL ) {
        factors_out.AddMatrices( this->A, this->nMatrices, device );
    }

    factors_out.FactorizeMIC();
}

void SparseSolverMIC::Solve( SEQ_VECTOR <double> ** rhs_sol) {

    if (!m_factorized) {
        //std::cout << "NOT INITIALIED\n\n"<<std::endl;
        std::stringstream ss;
        ss << "Solve -> rank: ";// << esinfo::mpi::MPIrank; // MPIrank link problem
        Factorization(ss.str());
    }
    esint offset = 0;
    
    for (esint i = 0; i < nMatrices; i++) {
        memcpy(vectors + offset, &(rhs_sol[i]->at(0)), m_Kplus_size[i] * sizeof(double));
        tmp_sol_d1[i] = vectors + offset;
        offset += m_Kplus_size[i];
            }
#pragma offload_transfer target(mic:device) in(vectors : length(offset) alloc_if(0) free_if(0)) \
    in(tmp_sol_d1 : length(nMatrices) alloc_if(0) free_if(0))
    if( USE_FLOAT ) {
        for (esint i = 0; i < nMatrices; i++) {
            for (esint j = 0; j < m_Kplus_size[i]; j++)
                tmp_sol_fl1[i][j] = (float)(rhs_sol[i])->at(j);
        }
    }
    double ** tmpVecPointers = new double*[nMatrices];

    if ( USE_FLOAT ) {
        for (esint i = 0 ; i < nMatrices; i++) {
            float *tmpPointer = tmp_sol_fl1[i];
            esint tmpLength = m_Kplus_size[i];
#pragma offload_transfer target(mic:device) \
            in(tmpPointer : length(tmpLength) alloc_if(0) free_if(0))
            //in(tmp_sol_fl1 : length(0) alloc_if(0) free_if(0))
            //{
            //    tmp_sol_fl1[i] = tmpPoesinter;
            //}
        }
    } else {
        for (esint i = 0; i < nMatrices; i++) {
            tmpVecPointers[i] = &(rhs_sol[i]->at(0));
        }
        //#pragma offload_transfer taget(mic:device) \
        //in(tmpVecPoesinters : length(nMatrices) alloc_if(1) free_if(0))
        for (esint i = 0; i < nMatrices;i++) {
            double * tmpPointer2 = tmpVecPointers[i];
            esint tmpLength = m_Kplus_size[i];
            double * buffer = tmp_sol_d1[i];
//#pragma offload_transfer target(mic:device) \
//          :  in(tmpPointer2 : length(tmpLength) alloc_if(0) free_if(0) into(buffer))
            //{
            //    tmpVecPoesinters[i] = tmpPoesinter2;
            //}

        }
    }

#pragma offload target(mic:device) \
    in(rows : length(0) alloc_if(0) free_if(0) ) \
    in(cols : length(0) alloc_if(0) free_if(0) ) \
    in(nnz : length(0) alloc_if(0) free_if(0) ) \
    in(CSR_I_row_indices_size : length(0) alloc_if(0) free_if(0)) \
    in(CSR_J_col_indices_size : length(0) alloc_if(0) free_if(0)) \
    in(CSR_I_row_indices : length(0) alloc_if(0) free_if(0)) \
    in(CSR_J_col_indices : length(0) alloc_if(0) free_if(0)) \
    in(CSR_V_values :length(0) alloc_if(0) free_if(0)) \
    in(pt : length(0) alloc_if(0) free_if(0)) \
    in(iparm : length(0) alloc_if(0) free_if(0)) \
    in(dparm : length(0) alloc_if(0) free_if(0)) \
    in(perm : length(0) alloc_if(0) free_if(0)) \
    in(error : length(0) alloc_if(0) free_if(0)) \
    in(m_Kplus_size : length(0) alloc_if(0) free_if(0)) \
    in(tmp_sol_d1 : length(0) alloc_if(0) free_if(0)) \
    in(tmp_sol_d2 : length(0) alloc_if(0) free_if(0)) \
    in(vectors : length(0) alloc_if(0) free_if(0)) \ 
    in(vectors_out : length(0) alloc_if(0) free_if(0)) \
    in(this : length(0) alloc_if(0) free_if(0))
    {
        esint offset = 0;
for (esint i = 0; i < nMatrices; i++){
    tmp_sol_d1[i] = &vectors[offset];
    tmp_sol_d2[i] = &vectors_out[offset];
    offset+=m_Kplus_size[i];
}

        double ddum   = 0;			/* Double dummy */
        MKL_INT idum  = 0;			/* Integer dummy. */
        MKL_INT n_rhs = 1;
        msglvl        = 0;
        mnum = 1;
        maxfct = 1;
        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */

        //iparm[24] = 1;		// Parallel forward/backward solve control. - 1 - Intel MKL PARDISO uses the sequential forward and backward solve.

#pragma omp parallel //num_threads(21)
        {
            esint myPhase = 331;
#pragma omp for schedule(dynamic)
            for (esint i = 0; i < nMatrices; i++) {
                myPhase=331;
                MKL_INT ip5backup = iparm[i][5];
                iparm[i][5] = 0;
                if (USE_FLOAT) {
                    myPhase=33;
                    pardiso (pt[i], &maxfct, &mnum, &mtype, &myPhase,
                            &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i], &idum, &n_rhs, iparm[i], &msglvl, &tmp_sol_fl1[i], &tmp_sol_fl2[i], &error[i]);


                } else {
myPhase = 33;
                    PARDISO (pt[i], &maxfct, &mnum, &mtype, &myPhase,
                            &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i], NULL, &n_rhs, iparm[i], &msglvl, tmp_sol_d1[i], tmp_sol_d2[i], &error[i]);



                    //myPhase = 332;
                    //PARDISO (pt[i], &maxfct, &mnum, &mtype, &myPhase,
                    //        &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i], &idum, &n_rhs, iparm[i], &msglvl, tmp_sol_d1[i], tmp_sol_d2[i], &error[i]);

                }
                iparm[i][5] = ip5backup;
            }
        }
        bool err = false;
        for (esint i = 0; i < nMatrices; i++) {
            if (error[i]!=0) {
                err = true;
                ////ESINFO(ERROR) << "ERROR during solution: " << error[i] << ", matrix " << i;
            }
        }
        if (err)
        {
            exit (3);
        }

        if (!keep_factors) {
            /* -------------------------------------------------------------------- */
            /* .. Termination and release of memory. */
            /* -------------------------------------------------------------------- */
            /*
               phase = -1;	 		
               MKL_INT nRhs = 1;
               for (esint i =0 ; i < nMatrices; i++ ) {
               pardiso (pt[i], &maxfct, &mnum, &mtype, &phase,
               &rows[i], &ddum, CSR_I_row_indices[i], CSR_J_col_indices[i], &idum, &nRhs,
               iparm[i], &msglvl, &ddum, &ddum, &error[i]);
               }
               initialized = false;
               */
        }


    }

    // send data from buffers on mic back to CPU user provided arrays
    if (USE_FLOAT) {
        for (esint i = 0; i < nMatrices; i++) {
            float * tmp = tmp_sol_fl2[i];
#pragma offload_transfer target(mic:device) \
            out(tmp : length(m_Kplus_size[i]) alloc_if(0) free_if(0))
            for (esint j = 0; j < m_Kplus_size[i]; j++) {
                rhs_sol[i]->at(j) = (double) tmp[j];
            }
        }
    } else {
#pragma offload_transfer target(mic:device) \
        out(vectors_out : length(offset) alloc_if(0) free_if(0))
        offset = 0;
        for (esint i = 0; i < nMatrices; i++) {
            memcpy(&(rhs_sol[i]->at(0)), vectors_out + offset, m_Kplus_size[i] * sizeof(double));
            offset += m_Kplus_size[i];
   //double * tmp = tmp_sol_d2[i];
            //double * output = tmpVecPointers[i];
//#pragma offload_transfer target(mic:device) \
//            out(tmp : length(m_Kplus_size[i]) alloc_if(0) free_if(0) into(output))
        }
    }
}

void SparseSolverMIC::Create_SC(
        DenseMatrixPack & SC_out,
        MKL_INT *SC_sizes,
        esint generate_symmetric_sc_1_generate_general_sc_0
        ) {
    ////ESINFO(PROGRESS3) << "Creating Schur complements";


    // data to be transfered to MIC
    esint * SC_out_rows = SC_out.rows;
    esint * SC_out_cols = SC_out.cols;
    long * SC_out_lengths = SC_out.lengths;
    long * SC_out_offsets = SC_out.offsets;
    long * SC_out_rowOffsets = SC_out.rowOffsets;
    long * SC_out_colOffsets = SC_out.colOffsets;
    long SC_out_totalCols = SC_out.totalCols;
    long SC_out_totalRows = SC_out.totalRows;
    double * SC_out_matrices = SC_out.matrices;
    esint matricesSize = SC_out.preallocSize - SC_out.freeSpace;
    bool * SC_out_packed = SC_out.packed;
    int setMsglvl = Info::report(LIBRARIES) ? 1 : 0;

#pragma offload target(mic : device ) \
    in ( rows : length(0) alloc_if(0) free_if(0) ) \
    in( SC_sizes  : length(nMatrices) alloc_if(1) free_if(1) ) \
    in( nnz : length(0) alloc_if(0) free_if(0))  \
    in( perm : length(0) alloc_if(0) free_if(0)) \
    in( CSR_V_values : length(0) alloc_if(0) free_if(0)  ) \
    in( CSR_I_row_indices  : length(0) alloc_if(0) free_if(0) ) \
    in( CSR_J_col_indices : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(SC_out_cols : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(SC_out_offsets : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(SC_out_matrices : length(matricesSize) alloc_if(1) free_if(0)) \
    in(SC_out_packed : length(nMatrices) alloc_if(1) free_if(0)) \
    in(SC_out_rowOffsets : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(SC_out_colOffsets : length(nMatrices) alloc_if(1) free_if (0)) \
    in(SC_out_lengths : length(nMatrices) alloc_if(1) free_if(0)) if(1) \
    in(this : length(0) alloc_if(0) free_if(0))
    {

        long maxSize = 0;
        for (esint i = 0; i < nMatrices; i++) {
            if (SC_sizes[i]*SC_sizes[i] > maxSize) {
                maxSize = SC_sizes[i]*SC_sizes[i];
            }
        }

        double ** matrixPerThread = new double*[nMatrices];
#pragma omp parallel
        {
            esint myRank = omp_get_thread_num();
            matrixPerThread[myRank] = (double*) _mm_malloc(maxSize * sizeof(double), 64);
        }
        ////ESINFO(PROGRESS3) << "start";
#pragma omp parallel
        {
            esint myRank = omp_get_thread_num();


#pragma omp for //schedule(static,1)
            for (esint i = 0 ; i < nMatrices; i++) {
                for (esint j = 0; j < 64; j++) {
                    iparm[i][j] = 0;
                }

                esint 	mtype = 2;
                double ddum;

                iparm[i][2] = 0;

                iparm[i][1-1] = 1;         /* No solver default */
                iparm[i][2-1] = 2;         /* Fill-in reordering from METIS */
                iparm[i][10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
                iparm[i][11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
                iparm[i][13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
                iparm[i][14-1] = 0;        /* Output: Number of perturbed pivots */
                iparm[i][18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
                iparm[i][19-1] = -1;       /* Output: Mflops for LU factorization */
                iparm[i][36-1] = 1;        /* Use Schur complement */

                maxfct = 1;           /* Maximum number of numerical factorizations. */
                mnum = 1;             /* Which factorization to use. */
                msglvl = setMsglvl;          /* Print statistical information in file */
                error = 0;            /* Initialize error flag */


                // for reordering and factorization
                for (esint j = 0; j < rows[i] - SC_sizes[i]; j++)
                    perm[i][j] = 0;

                for (esint j = rows[i] - SC_sizes[i]; j < rows[i]; j++)
                    perm[i][j] = 1;
                m_nRhs = 0;

                /* -------------------------------------------------------------------- */
                /* .. Numerical factorization. */
                /* -------------------------------------------------------------------- */

                double * SC_out_data;
                if (SC_out_packed[i]) {
                    SC_out_data = matrixPerThread[myRank];//(double*) malloc(K_b_tmp_rows[i] * K_b_tmp_rows[i]*sizeof(double));
                } else {
                    SC_out_data = SC_out_matrices + SC_out_offsets[i];
                }
                phase = 12;
                PARDISO (pt[i], &maxfct, &mnum, &mtype, &phase,
                        &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i],
                        perm[i], &m_nRhs,
                        iparm[i], &msglvl, &ddum, SC_out_data, &error[i]);


                for (esint j = 0; j < SC_sizes[i]*SC_sizes[i]; j++)
                    SC_out_data[j] = (-1.0)*SC_out_data[j];

                if (SC_out_packed[i]) {
                    for (long j = 0; j < SC_out_cols[ i ]; j++ ) {
                        long pos = ((double)( 1.0 +j ))*((double) j)/2.0;
                        memcpy( SC_out_matrices + SC_out_offsets[ i ] + pos,
                                SC_out_data + j * SC_out_rows[ i ], sizeof( double ) * (j+1) );
                    }
                }
                /* -------------------------------------------------------------------- */
                /* .. Termination and release of memory. */
                /* -------------------------------------------------------------------- */
                phase = -1;           /* Release internal memory. */

                PARDISO (pt[i], &maxfct, &mnum, &mtype, &phase,
                        &rows[i], CSR_V_values[i], CSR_I_row_indices[i], CSR_J_col_indices[i],
                        perm[i], &m_nRhs,
                        iparm[i], &msglvl, &ddum, SC_out_data, &error[i]);

                initialized = false;
            }
        }
        bool err=false;
        for (esint i = 0; i < nMatrices; i++) {
            if (error[i]!=0) {
                err=true;
            }
        }
        if ( err )
        {
            ////ESINFO(ERROR) << "ERROR during numerical factorization";
            exit (2);
        } else {
            initialized = true;
        }


        ////ESINFO(PROGRESS3) << "end";
    }
}

void SparseSolverMIC::Create_SC_w_Mat(
        SparseMatrix** K_in,
        SparseMatrix** B_in,
        DenseMatrixPack & SC_out,
        esint nMatrices,
        esint generate_symmetric_sc_1_generate_general_sc_0,
        esint device ) {
    if (!SC_out.areDataOnMIC()) {
        SC_out.CopyToMIC();
    }
    // data to be transfered to MIC
    esint * K_in_rows = (esint*) _mm_malloc(nMatrices * sizeof(esint), 64); //new esint[nMatrices];
    esint * K_b_tmp_rows = (esint*) _mm_malloc(nMatrices * sizeof(esint), 64);//new esint[nMatrices];
    esint * K_sc1_rows = (esint*) _mm_malloc(nMatrices * sizeof(esint), 64);//new esint[nMatrices];
    esint * nnz = (esint*) _mm_malloc(nMatrices * sizeof(esint), 64);//new esint[nMatrices];

    double ** K_sc1_CSR_V_values = (double**) _mm_malloc(nMatrices * sizeof(double*), 64);// new double*[nMatrices];
    esint ** K_sc1_CSR_I_row_indices = (esint**) _mm_malloc(nMatrices * sizeof(esint*), 64); //new esint*[nMatrices];
    esint ** K_sc1_CSR_J_col_indices = (esint**) _mm_malloc(nMatrices * sizeof(esint*), 64);//new esint*[nMatrices];
    //double ** SC_out_dense_values = new double*[nMatrices];
    esint * SC_out_rows = SC_out.rows;
    esint * SC_out_cols = SC_out.cols;
    long * SC_out_lengths = SC_out.lengths;
    long * SC_out_offsets = SC_out.offsets;
    long * SC_out_rowOffsets = SC_out.rowOffsets;
    long * SC_out_colOffsets = SC_out.colOffsets;
    double * SC_out_mic_x_in = SC_out.mic_x_in;
    double * SC_out_mic_y_out = SC_out.mic_y_out;
    long SC_out_totalCols = SC_out.totalCols;
    long SC_out_totalRows = SC_out.totalRows;
    double * SC_out_matrices = SC_out.matrices;
    esint matricesSize = SC_out.preallocSize - SC_out.freeSpace;
    bool * SC_out_packed = SC_out.packed;

    // find the biggest SC matrix and preallocate output array for PARDISO
    for ( esint i = 0 ; i < nMatrices; ++i ) {


        // assemble the whole matrix [K, B; Bt, eye]
        SparseMatrix K_sc1;
        SparseMatrix Sc_eye;
        SparseMatrix K_b_tmp = *(B_in[i]);

        K_b_tmp.MatTranspose();
        Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);
        K_sc1 = *(K_in[i]);
        K_sc1.MatTranspose();
        K_sc1.MatAppend(K_b_tmp);
        K_sc1.MatTranspose();
        K_sc1.MatAppend(Sc_eye);

        K_in_rows[ i ] = K_in[i]->rows;
        K_b_tmp_rows[i] = K_b_tmp.rows;
        K_sc1_rows[i] = K_sc1.rows;
        nnz[i] = K_sc1.nnz;

        K_sc1_CSR_V_values[i] = (double*) _mm_malloc(nnz[i] * sizeof(double), 64); // new double[ nnz[i] ];
        memcpy(K_sc1_CSR_V_values[i], &K_sc1.CSR_V_values[0], nnz[i] * sizeof(double) );
        K_sc1_CSR_I_row_indices[i] = (esint*) _mm_malloc( (K_sc1_rows[i] +1) * sizeof(esint), 64);// new esint[K_sc1_rows[i] + 1];
        memcpy(K_sc1_CSR_I_row_indices[i], &K_sc1.CSR_I_row_indices[0], ( K_sc1_rows[i] + 1) * sizeof(esint));
        //K_sc1_CSR_I_row_indices[i][K_sc1_rows[i]] = nnz[i];
        K_sc1_CSR_J_col_indices[i] = (esint*) _mm_malloc(nnz[i] * sizeof(esint), 64); //new esint[nnz[i]];
        memcpy(K_sc1_CSR_J_col_indices[i], &K_sc1.CSR_J_col_indices[0], nnz[i] * sizeof(esint));

    }
    // transfer data
#pragma offload_transfer target(mic : device ) \
    in( K_in_rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in( K_b_tmp_rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in( K_sc1_rows : length(nMatrices) alloc_if(1) free_if(0) ) \
    in( nnz : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(K_sc1_CSR_V_values : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(K_sc1_CSR_I_row_indices : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(K_sc1_CSR_J_col_indices : length(nMatrices) alloc_if(1) free_if(0) ) \
    in(this : alloc_if(1) free_if(0) )
    // #pragma omp parallel for num_threads(24)
    for (esint i = 0 ; i < nMatrices; i++) {
        double * valPointer = K_sc1_CSR_V_values[i];
        esint * iPointer = K_sc1_CSR_I_row_indices[i];
        esint * jPointer = K_sc1_CSR_J_col_indices[i];
#pragma offload target(mic:device) \
        in( valPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( iPointer : length(K_sc1_rows[i] + 1) alloc_if(1) free_if(0)) \
        in( jPointer : length(nnz[i]) alloc_if(1) free_if(0)) \
        in( K_sc1_CSR_V_values : length(0) alloc_if(0) free_if(0) ) \
        in( K_sc1_CSR_I_row_indices : length(0) alloc_if(0) free_if(0) ) \
        in( K_sc1_CSR_J_col_indices : length(0) alloc_if(0) free_if(0) ) \
        in(this : length(0) alloc_if(0) free_if(0))
        {
            K_sc1_CSR_V_values[i] = valPointer;
            K_sc1_CSR_I_row_indices[i] = iPointer;
            K_sc1_CSR_J_col_indices[i] = jPointer;
        }
    }
#pragma offload target(mic : device ) \
    in ( K_in_rows : length(0) alloc_if(0) free_if(1) ) \
    in( K_b_tmp_rows : length(0) alloc_if(0) free_if(1) ) \
    in( K_sc1_rows  : length(0) alloc_if(0) free_if(1) ) \
    in( nnz : length(0) alloc_if(0) free_if(1))  \
    in( K_sc1_CSR_V_values : length(0) alloc_if(0) free_if(0)  ) \
    in(K_sc1_CSR_I_row_indices  : length(0) alloc_if(0) free_if(0) ) \
    in( K_sc1_CSR_J_col_indices : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_rows : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_cols : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_offsets : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_matrices : length(0) alloc_if(0) free_if(0)) \
    in(SC_out_packed : length(0) alloc_if(0) free_if(0)) \
    in(SC_out_rowOffsets : length(0) alloc_if(0) free_if(0) ) \
    in(SC_out_colOffsets : length(0) alloc_if(0) free_if (0)) \
    in(SC_out_mic_x_in : length(0) alloc_if(0) free_if(0)) \
    in(SC_out_mic_y_out : length(0) alloc_if(0) free_if(0)) \
    in(SC_out_lengths : length(0) alloc_if(0) free_if(0))  in(this : length(0) alloc_if(0) free_if(0)) if(1)
    {
        long maxSize = 0;
        for (esint i = 0; i < nMatrices; i++) {
            if (K_b_tmp_rows[i]*K_b_tmp_rows[i] > maxSize) {
                maxSize = K_b_tmp_rows[i]*K_b_tmp_rows[i];
            }
        }
        esint nThreads = 1;
#pragma omp parallel
        {
#pragma omp single
           nThreads = omp_get_num_threads();
            }
        double ** matrixPerThread = new double*[nThreads];
#pragma omp parallel
        {
            esint myRank = omp_get_thread_num();
            // to overcome competition among threads allocate buffer one by one
#pragma omp critical
            matrixPerThread[myRank] = (double*) _mm_malloc(maxSize * sizeof(double), 64);

            /* Internal solver memory pointer pt, */
            /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
            /* or void *pt[64] should be OK on both architectures */
            void *pt[64];

            /* Pardiso control parameters. */
            esint 	iparm[64];
            double  dparm[65];
            esint 	maxfct, mnum, phase, error;

            /* Auxiliary variables. */
            esint 	j;
            double 	ddum;			/* Double dummy */
            esint 	idum;			/* Integer dummy. */
            esint 	solver;

            /* -------------------------------------------------------------------- */
            /* .. Setup Pardiso control parameters. */
            /* -------------------------------------------------------------------- */
            for (j = 0; j < 64; j++) {
                iparm[j] = 0;
            }

            /* -------------------------------------------------------------------- */
            /* .. Initialize the esinternal solver memory poesinter. This is only */
            /* necessary for the FIRST call of the PARDISO solver. */
            /* -------------------------------------------------------------------- */
            for (j = 0; j < 64; j++)
                pt[j] = 0;

            esint 	mtype = 2;

            iparm[2] = 0;

            iparm[1-1] = 1;         /* No solver default */
            iparm[2-1] = 2;         /* Fill-in reordering from METIS */
            iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
            iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
            iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
            iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
            iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
            iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
            iparm[36-1] = 1;        /* Use Schur complement */

            maxfct = 1;           /* Maximum number of numerical factorizations. */
            mnum = 1;             /* Which factorization to use. */
            msglvl = 0;          /* Print statistical information in file */
            error = 0;            /* Initialize error flag */

#pragma omp for //schedule(static,1)
            for (esint i = 0 ; i < nMatrices; i++) {
                // for reordering and factorization
                esint * perm = new esint[K_sc1_rows[i]];
                for (esint j = 0; j < K_in_rows[i]; j++)
                    perm[j] = 0;

                for (esint j = K_in_rows[i]; j < K_sc1_rows[i]; j++)
                    perm[j] = 1;
                esint nrhs = 0;

                /* -------------------------------------------------------------------- */
                /* .. Numerical factorization. */
                /* -------------------------------------------------------------------- */

                double * SC_out_data;
                if (SC_out_packed[i]) {
                    SC_out_data = matrixPerThread[myRank];//(double*) malloc(K_b_tmp_rows[i] * K_b_tmp_rows[i]*sizeof(double));
                } else {
                    SC_out_data = SC_out_matrices + SC_out_offsets[i];
                }

                phase = 12;

                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                        &K_sc1_rows[i],
                        K_sc1_CSR_V_values[i], K_sc1_CSR_I_row_indices[i], K_sc1_CSR_J_col_indices[i],
                        &perm[0], &nrhs,
                        iparm, &msglvl, &ddum, SC_out_data, &error);

                for (esint j = 0; j < K_b_tmp_rows[i]*K_b_tmp_rows[i]; j++)
                    SC_out_data[j] = (-1.0)*SC_out_data[j];

                if (SC_out_packed[i]) {
                    for (long j = 0; j < SC_out_cols[ i ]; j++ ) {
                        long pos = ((double)( 1.0 +j ))*((double) j)/2.0;
                        memcpy( SC_out_matrices + SC_out_offsets[ i ] + pos,
                                SC_out_data + j * SC_out_rows[ i ], sizeof( double ) * (j+1) );
                    }
                }
                if ( error != 0 )
                {
                    ////ESINFO(ERROR) << "ERROR during numerical factorization: " << error;
                    exit (2);
                } else {
                    initialized = true;
                }

                /* -------------------------------------------------------------------- */
                /* .. Termination and release of memory. */
                /* -------------------------------------------------------------------- */
                phase = -1;           /* Release internal memory. */


                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                        &K_sc1_rows[i],
                        K_sc1_CSR_V_values[i], K_sc1_CSR_I_row_indices[i], K_sc1_CSR_J_col_indices[i],
                        &perm[0], &nrhs,
                        iparm, &msglvl, &ddum, SC_out_data, &error);


                initialized = false;

                delete [] perm;
            }
            _mm_free(matrixPerThread[myRank]);
        }
    }
        // remember to clear also MICs memory
    for (esint i = 0 ; i < nMatrices; i++) {
        double * valPointer = K_sc1_CSR_V_values[i];
        esint * iPointer = K_sc1_CSR_I_row_indices[i];
        esint * jPointer = K_sc1_CSR_J_col_indices[i];
#pragma offload_transfer target(mic:device) \
        nocopy( valPointer : alloc_if(0) free_if(1)) \
        nocopy( iPointer :  alloc_if(0) free_if(1)) \
        nocopy( jPointer :  alloc_if(0) free_if(1)) 
    }
  #pragma offload_transfer target(mic : device ) \
    nocopy(K_sc1_CSR_V_values :alloc_if(0) free_if(1) ) \
    nocopy(K_sc1_CSR_I_row_indices :  alloc_if(0) free_if(1) ) \
    nocopy(K_sc1_CSR_J_col_indices : alloc_if(0) free_if(1) )

 for ( esint i = 0 ; i < nMatrices; ++i ) {
        _mm_free(K_sc1_CSR_V_values[i]);
        _mm_free(K_sc1_CSR_I_row_indices[i]);
        _mm_free(K_sc1_CSR_J_col_indices[i]);
    }
    _mm_free(K_sc1_CSR_V_values);
    _mm_free(K_sc1_CSR_I_row_indices);
    _mm_free(K_sc1_CSR_J_col_indices);
    _mm_free(K_in_rows);
    _mm_free(K_b_tmp_rows);
    _mm_free(K_sc1_rows);
    _mm_free(nnz);
    

}
