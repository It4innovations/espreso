
#ifdef HAVE_SUPERLU_DIST

#include "wrappers/superlu_dist/operations/solver_csx.superlu_dist.h"

#include "wrappers/superlu_dist/superlu_dist_common.h"

#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_csx.h"
#include "math/operations/complete_csx_csy_map.h"
#include "math/operations/convert_csx_csy_map.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/convert_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



namespace copied_code {
    void my_superlu_gridmap(
                MPI_Comm Bcomm, /* The base communicator upon which
                        the new grid is formed. */
                int nprow,
                int npcol,
                int usermap[], /* usermap(i,j) holds the process
                        number to be placed in {i,j} of
                        the process grid.  */
                int ldumap,    /* The leading dimension of the
                        2D array usermap[].  */
                gridinfo_t *grid)
    {
        MPI_Group mpi_base_group, superlu_grp;
        int Np = nprow * npcol, mycol, myrow;
        int *pranks;
        int i, j, info;

    #if 0 // older MPI doesn't support complex in C    
        /* Create datatype in C for MPI complex. */
        if ( SuperLU_MPI_DOUBLE_COMPLEX == MPI_DATATYPE_NULL ) {
            MPI_Type_contiguous( 2, MPI_DOUBLE, &SuperLU_MPI_DOUBLE_COMPLEX );
            MPI_Type_commit( &SuperLU_MPI_DOUBLE_COMPLEX );
        }
    #endif
        
        /* Check MPI environment initialization. */
        MPI_Initialized( &info );
        if ( !info )
            ABORT("C main program must explicitly call MPI_Init()");

        grid->nprow = nprow;
        grid->npcol = npcol;

        /* Make a list of the processes in the new communicator. */
        pranks = (int *) SUPERLU_MALLOC(Np*sizeof(int));
        for (j = 0; j < npcol; ++j)
        for (i = 0; i < nprow; ++i)
            pranks[i*npcol+j] = usermap[j*ldumap+i];
        
        /*
        * Form MPI communicator for all.
        */
        /* Get the group underlying Bcomm. */
        MPI_Comm_group( Bcomm, &mpi_base_group );
        /* Create the new group. */
        MPI_Group_incl( mpi_base_group, Np, pranks, &superlu_grp );
        /* Create the new communicator. */
        /* NOTE: The call is to be executed by all processes in Bcomm,
        even if they do not belong in the new group -- superlu_grp.
        The function returns MPI_COMM_NULL to processes that are not in superlu_grp. */
        MPI_Comm_create( Bcomm, superlu_grp, &grid->comm );

        /* Bail out if I am not in the group "superlu_grp". */
        if ( grid->comm == MPI_COMM_NULL ) {
            // grid->comm = Bcomm;  do not need to reassign to a valid communicator
            grid->iam = -1;
            //SUPERLU_FREE(pranks);
            //return;
            goto gridmap_out;
        }

        MPI_Comm_rank( grid->comm, &(grid->iam) );
        myrow = grid->iam / npcol;
        mycol = grid->iam % npcol;

        /*
        * Form MPI communicator for myrow, scope = COMM_ROW.
        */
    #if 0
        for (i = 0; i < npcol; ++i) pranks[i] = myrow*npcol + i;
        MPI_Comm_group( grid->comm, &superlu_grp );          /* Find all's group */
        MPI_Group_incl( superlu_grp, npcol, pranks, &grp );  /* Form new group */
        MPI_Comm_create( grid->comm, grp, &grid->rscp.comm );/* Create new comm */
    #else
        MPI_Comm_split(grid->comm, myrow, mycol, &(grid->rscp.comm));
    #endif

        /*
        * Form MPI communicator for mycol, scope = COMM_COLUMN.
        */
    #if 0
        for (i = 0; i < nprow; ++i) pranks[i] = i*npcol + mycol;
        MPI_Group_incl( superlu_grp, nprow, pranks, &grp );  /* Form new group */
        MPI_Comm_create( grid->comm, grp, &grid->cscp.comm );/* Create new comm */
    #else
        MPI_Comm_split(grid->comm, mycol, myrow, &(grid->cscp.comm));
    #endif

        grid->rscp.Np = npcol;
        grid->rscp.Iam = mycol;
        grid->cscp.Np = nprow;
        grid->cscp.Iam = myrow;

    #if 0
        {
        int tag_ub;
        if ( !grid->iam ) {
            MPI_Comm_get_attr(Bcomm, MPI_TAG_UB, &tag_ub, &info);
            printf("MPI_TAG_UB %d\n", tag_ub);
            /* returns 4295677672
            In reality it is restricted to no greater than 16384. */
        }
        exit(0);
        }
    #endif

    gridmap_out:        
        SUPERLU_FREE(pranks);
        MPI_Group_free(&superlu_grp);
        MPI_Group_free(&mpi_base_group);
        
    } /* superlu_gridmap */

    void my_superlu_gridinit(MPI_Comm Bcomm, /* The base communicator upon which
                        the new grid is formed. */
                int nprow, int npcol, gridinfo_t *grid)
    {
        int Np = nprow * npcol;
        int *usermap;
        int i, j, info;

        /* Make a list of the processes in the new communicator. */
        usermap = (int*)SUPERLU_MALLOC(Np*sizeof(int));
        for (j = 0; j < npcol; ++j)
        for (i = 0; i < nprow; ++i) usermap[j*nprow+i] = i*npcol+j;
        
        /* Check MPI environment initialization. */
        MPI_Initialized( &info );
        if ( !info )
            ABORT("C main program must explicitly call MPI_Init()");

        MPI_Comm_size( Bcomm, &info );
        if ( info < Np ) {
            printf("Number of processes %d is smaller than NPROW * NPCOL %d", info, Np);
            exit(-1);
        }

        copied_code::my_superlu_gridmap(Bcomm, nprow, npcol, usermap, nprow, grid);
        
        SUPERLU_FREE(usermap);
        
    // #ifdef GPU_ACC
        /* Binding each MPI to a GPU device */
        char *ttemp;
        ttemp = getenv ("SUPERLU_BIND_MPI_GPU");

        if (ttemp) {
            int devs, rank;
            MPI_Comm_rank(Bcomm, &rank); // MPI_COMM_WORLD??
            gpuGetDeviceCount(&devs);  // Returns the number of compute-capable devices
            int device_id = (int)(rank/get_mpi_process_per_gpu ()) % (devs); //YL: allow multiple MPIs per GPU
            gpuSetDevice(device_id); // Set device to be used for GPU executions

            int get_cur_dev;
            gpuGetDevice(&get_cur_dev);
            printf("** MPI rank %d, gpu=%d **\n",rank,get_cur_dev);
            fflush(stdout);
        }
    // #endif
    }

}



static int num_active_intances = 0;



template<typename T, typename I>
struct solver_csx_superlu_dist_data
{
    struct config {
        bool use_gpu = false;
    } cfg;
    MatrixCsxData_new<T,I> A_csc_full;
    MatrixCsxData_new<T,I> A_into_slu;
    complete_csx_csy_map<T,I> op_A_complete;
    convert_csx_csy_map<T,I> op_A_reorder;
    copy_csx<T,I> op_A_copy;
    gridinfo_t grid;
    NCformat slu_A_store;
    SuperMatrix slu_A;
    dScalePermstruct_t scale_permstruct;
    dLUstruct_t lu_struct;
    SuperLUStat_t stat;
    std::vector<double> berr;
    superlu_dist_options_t options;
    bool need_complete;
    bool need_reorder;
    bool called_factorize_numeric = false;
};



template<typename T, typename I>
solver_csx_superlu_dist<T,I>::solver_csx_superlu_dist()
{
    if(!std::is_same_v<T,double>) eslog::error("only T=double is supported\n");

    #pragma omp critical(solver_csx_superlu_dist)
    {
        num_active_intances++;
        if(num_active_intances > 1) {
            eslog::error("only one instance of solver_csx_superlu_dist can be active at once\n");
        }
    }

    data = std::make_unique<solver_csx_superlu_dist_data<T,I>>();
}



template<typename T, typename I>
solver_csx_superlu_dist<T,I>::~solver_csx_superlu_dist()
{
    dScalePermstructFree(&data->scale_permstruct);
    dLUstructFree(&data->lu_struct);
    PStatFree(&data->stat);
    superlu_gridexit(&data->grid);

    #pragma omp critical(solver_csx_superlu_dist)
    {
        num_active_intances--;
    }
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_factorize_symbolic()
{
    data->need_complete = is_uplo(A->prop.uplo);
    data->need_reorder = !data->need_complete && (A->order != 'C');

    if(data->need_complete || data->need_reorder) {
        data->A_csc_full.set(A->nrows, A->ncols, 2 * A->nnz - A->nrows, 'C', AllocatorCPU_new::get_singleton());
        data->A_csc_full.prop = A->prop;
        data->A_csc_full.prop.uplo = 'F';
        data->A_csc_full.alloc();
    }

    if(data->need_complete) {
        data->op_A_complete.set_matrix_src(A);
        data->op_A_complete.set_matrix_dst(&data->A_csc_full);
        data->op_A_complete.set_conj(true);
        data->op_A_complete.perform_pattern();
    }
    if(data->need_reorder) {
        data->op_A_reorder.set_matrix_src(A);
        data->op_A_reorder.set_matrix_dst(&data->A_csc_full);
        data->op_A_reorder.perform_pattern();
    }

    // because superlu modifies the data I give it, I need to create a copy
    if(data->need_complete || data->need_reorder) {
        data->A_into_slu.set(data->A_csc_full.nrows, data->A_csc_full.ncols, data->A_csc_full.nnz, 'C', AllocatorCPU_new::get_singleton());
        data->A_into_slu.prop = data->A_csc_full.prop;
        data->A_into_slu.prop.uplo = 'F';
        data->A_into_slu.alloc();
        data->op_A_copy.set_matrix_src(&data->A_csc_full);
    }
    else {
        data->A_into_slu.set(A->nrows, A->ncols, A->nnz, 'C', AllocatorCPU_new::get_singleton());
        data->A_into_slu.prop = A->prop;
        data->A_into_slu.prop.uplo = 'F';
        data->A_into_slu.alloc();
        data->op_A_copy.set_matrix_src(A);
    }
    data->op_A_copy.set_matrix_dst(&data->A_into_slu);

    // superlu_gridinit(MPI_COMM_SELF, 1, 1, &data->grid);
    copied_code::my_superlu_gridinit(MPI_COMM_SELF, 1, 1, &data->grid);

    data->slu_A_store.nnz = data->A_into_slu.nnz;
    data->slu_A_store.nzval = data->A_into_slu.vals;
    data->slu_A_store.rowind = data->A_into_slu.idxs;
    data->slu_A_store.colptr = data->A_into_slu.ptrs;

    data->slu_A.Stype = SLU_NC; // CSC
    data->slu_A.Dtype = SLU_D; // double
    data->slu_A.Mtype = SLU_GE; // general
    data->slu_A.nrow = A->nrows;
    data->slu_A.ncol = A->ncols;
    data->slu_A.Store = &data->slu_A_store;

    dScalePermstructInit(A->nrows, A->ncols, &data->scale_permstruct);

    dLUstructInit(A->nrows, &data->lu_struct);

    PStatInit(&data->stat);

    set_default_options_dist(&data->options);
    data->options.Equil = NO;
    data->options.ParSymbFact = NO;
    data->options.ColPerm = METIS_AT_PLUS_A;
    data->options.RowPerm = NOROWPERM;
    data->options.ReplaceTinyPivot = NO;
    data->options.IterRefine = NOREFINE;
    data->options.Trans = NOTRANS;
    data->options.SymPattern = YES;
    data->options.PrintStat = NO;
    data->options.superlu_acc_offload = (data->cfg.use_gpu ? 1 : 0);
    data->options.superlu_num_gpu_streams = 1;

    // symbolic factorization is done lazily in numeric factorization
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_factorize_numeric()
{
    if(data->need_complete) {
        data->op_A_complete.perform_values();
    }
    if(data->need_reorder) {
        data->op_A_reorder.perform_values();
    }
    data->op_A_copy.perform();

    data->options.Fact = (data->called_factorize_numeric ? SamePattern_SameRowPerm : DOFACT);

    int info;
printf("GGG %p %d\n", &data->grid, data->grid.npcol);
// now it crashes here
    pdgssvx_ABglobal(&data->options, &data->slu_A, &data->scale_permstruct, nullptr, A->nrows, 0, &data->grid, &data->lu_struct, nullptr, &data->stat, &info);
    if(info != 0) eslog::error("superlu error, info = %d\n", info);

    data->called_factorize_numeric = true;
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    eslog::error("no support for permutation extraction\n");
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    eslog::error("no support for factor extraction\n");
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MatrixDenseView_new<T> rhs_mat;
    rhs_mat.set_view(rhs.size, 1, rhs.size, 'C', rhs.vals, rhs.ator);

    MatrixDenseView_new<T> sol_mat;
    sol_mat.set_view(sol.size, 1, sol.size, 'C', sol.vals, sol.ator);

    if(&rhs == &sol) {
        this->internal_solve(sol_mat, sol_mat);
    }
    else {
        this->internal_solve(rhs_mat, sol_mat);
    }
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order == 'R') {
        MatrixDenseData_new<T> tmp;
        tmp.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
        tmp.alloc();
        convert_dnx_dny<T>::do_all(&rhs, &tmp, false);
        this->internal_solve(tmp, tmp);
        convert_dnx_dny<T>::do_all(&tmp, &sol, false);
        return;
    }

    if(&rhs != &sol) {
        copy_dnx<T>::do_all(&rhs, &sol, false);
    }

    data->berr.resize(sol.ncols);

    data->options.Fact = FACTORED;

    int info;
    pdgssvx_ABglobal(&data->options, &data->slu_A, &data->scale_permstruct, nullptr, A->nrows, 0, &data->grid, &data->lu_struct, data->berr.data(), &data->stat, &info);
    if(info != 0) eslog::error("superlu error, info = %d\n", info);
}



template<typename T, typename I>
void solver_csx_superlu_dist<T,I>::internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order == 'R') {
        MatrixDenseData_new<T> tmp;
        tmp.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
        tmp.alloc();
        convert_csx_dny<T,I>::do_all(&rhs, &tmp);
        this->internal_solve(tmp, tmp);
        convert_dnx_dny<T>::do_all(&tmp, &sol, false);
    }
    else {
        convert_csx_dny<T,I>::do_all(&rhs, &sol);
        this->internal_solve(sol, sol);
    }
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_superlu_dist<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}

#endif
