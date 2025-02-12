
#include "w.wsmp.systemsolver.h"
#include "wrappers/mpi/communication.h"
#include "config/ecf/linearsolver/wsmp.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "physics/system/wsmpsystem.h"

#include <vector>
#include <cstring>

#ifdef HAVE_WSMP
extern "C" void pwgsmp_(int *n, int ia[], int ja[],
    double avals[], double b[], int *ldb,
    int *nrhs, double rmisc[], int iparm[], double dparm[]);

extern "C" void pwssmp_(int *ni, int iai[], int jai[],
    double avalsi[], double diagi[], int perm[], int invp[],
    double bi[], int *ldbi, int *nrhs, double *auxi, int *nauxi,
    int mrpi[], int iparm[], double dparm[]);

namespace espreso {
struct WSMPDataHolder {
    int m, n;
    int ldb, nrhs;
    double *solution;
    double *rhsValues;
    std::vector<esint> rowData;
    esint *colIndices;
    double *values;
    std::vector<esint> colData;
    std::vector<double> valData;
    int nvalues;
    double *rmisc;
    double *diag;
    int *mrp;
    int *perm;
    int *invp;
    int iparm[64];
    double dparm[64];
    int mtype;
    double aux;
    int naux;
};
}
#endif

using namespace espreso;

WSMPSystemSolver::WSMPSystemSolver(WSMPConfiguration &configuration, WSMPSolverData &data)
: configuration(configuration), _roffset(0), _nrows(0), _precision(1e-15), _data(data), _inner(NULL)
{
#ifndef HAVE_WSMP
    eslog::globalerror("ESPRESO run-time error: cannot call WSMP solver (the library with the solver is not linked).\n");
#endif
}

void WSMPSystemSolver::init()
{
#ifdef HAVE_WSMP
    _nrows = _data.K.nrows - _data.K.nhalo;
    _roffset = _nrows;

    _inner = new WSMPDataHolder();
    _inner->n = Communication::exscan(_roffset);
    _inner->nrhs = 1;
    _inner->m = _inner->n;
    _inner->ldb = _nrows;
    _inner->rmisc = new double[_inner->ldb * _inner->nrhs];
    _inner->diag = new double[_nrows];
    _inner->mrp = new int[_nrows];
    _inner->perm = new int[_inner->n];
    _inner->invp = new int[_inner->n];
    _inner->aux = 0;
    _inner->naux = 0;

    _inner->iparm[0] = 0;
#endif
}

void WSMPSystemSolver::update()
{
#ifdef HAVE_WSMP
    eslog::solver("     - ---- LINEAR SOLVER -------------------------------------------------------------- -\n");
    eslog::solver("     - | SOLVER ::     WSMP                      TYPE :: PARALLEL DIRECT SPARSE SOLVER | -\n");

    switch (_data.K.type) {
    case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
        eslog::solver("     - | MATRIX TYPE ::                               REAL SYMMETRIC POSITIVE DEFINITE | -\n");
        _inner->mtype = 2; break;
    case MatrixType::REAL_SYMMETRIC_INDEFINITE:
        eslog::solver("     - | MATRIX TYPE ::                                      REAL SYMMETRIC INDEFINITE | -\n");
        _inner->mtype = -2; break;
    case MatrixType::REAL_UNSYMMETRIC:
        eslog::solver("     - | MATRIX TYPE ::                      REAL NONSYMMETRIC (STRUCTURALY SYMMETRIC) | -\n");
        _inner->mtype = 1; break;
    }

    double start = eslog::time();

    const bool fillStructure = !_inner->rowData.size();

    // pick only upper triangle (since composer does not set correct
    // dirichlet in symmetric matrices)
    if (_data.K.type == MatrixType::REAL_SYMMETRIC_INDEFINITE || _data.K.type == MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
        _inner->iparm[1] = 3;
        if (fillStructure) {
            _inner->iparm[1] = 1;
            _inner->iparm[2] = 5;

            _inner->rowData.clear();
            _inner->colData.clear();
            _inner->rowData.reserve(_nrows + 1);
            _inner->rowData.push_back(1);
        }
        _inner->valData.clear();
        for (esint i = _data.K.nhalo; i < _data.K.nrows; i++) {
            for (esint c = _data.K.rows[i] - 1; c < _data.K.rows[i + 1] - 1; c++) {
                if (_roffset + i - _data.K.nhalo <= _data.K.cols[c] - 1) {
                    if (fillStructure) {
                        _inner->colData.push_back(_data.K.cols[c]);
                    }
                    _inner->valData.push_back(_data.K.vals[c]);
                }
            }
            if (fillStructure) {
                _inner->rowData.push_back(_inner->colData.size() + 1);
            }
        }
        _inner->colIndices = _inner->colData.data();
        _inner->values = _inner->valData.data();

        if (_data.K.type == MatrixType::REAL_SYMMETRIC_INDEFINITE)
        {
            // IPARM
            _inner->iparm[0] = 1;
            _inner->iparm[3] = 0;
            _inner->iparm[4] = _inner->iparm[5] = 1;
            _inner->iparm[6] = 3;
            for (int i = 7; i < 15; i++)
            { _inner->iparm[i] = 0; }
            // iparm[9], iparm[10] differ from default
            _inner->iparm[9] = 3;
            _inner->iparm[10] = 1;
            _inner->iparm[15] = _inner->iparm[17] = 1;
            _inner->iparm[16] = 0;
            _inner->iparm[18] = 0;
            _inner->iparm[19] = 0;
            for (int i = 24; i < 32; i++)
            { _inner->iparm[i] = 0; }
            // iparm[30] differs from default
            _inner->iparm[30] = 2;
            for (int i = 33; i < 63; i++)
            { _inner->iparm[i] = 0; }

            // DPARM
            _inner->dparm[5] = 2e-15;
            _inner->dparm[9] = 1e-18;
            _inner->dparm[10] = 0.001;
            _inner->dparm[11] = 4e-16;
            _inner->dparm[14] = 0;
            _inner->dparm[20] = 1e200;
            _inner->dparm[21] = 2e-8;
            _inner->dparm[30] = 2;
            for (int i = 36; i < 63; i++)
            { _inner->dparm[i] = 0; }
        }
    } else {
        _inner->iparm[1] = 2;

        /* Custom parameters
        // IPARM
        _inner->iparm[0] = 1;
        _inner->iparm[3] = 0;
        _inner->iparm[4] = 1;
        _inner->iparm[5] = _inner->iparm[6] = 3;
        _inner->iparm[7] = 0;
        _inner->iparm[8] = 0;
        _inner->iparm[9] = 1;
        _inner->iparm[10] = 1;
        _inner->iparm[11] = 0;
        _inner->iparm[14] = 25;
        _inner->iparm[15] = 1;
        for (int i = 16; i < 20; i++)
        { _inner->iparm[i] = 0; }
        _inner->iparm[20] = 1;
        _inner->iparm[24] = 0;
        for (int i = 26; i < 30; i++)
        { _inner->iparm[i] = 0; }
        _inner->iparm[30] = 1;
        _inner->iparm[31] = 0;
        _inner->iparm[33] = 10;
        for (int i = 34; i < 63; i++)
        { _inner->iparm[i] = 0; }

        // DPARM
        _inner->dparm[5] = 2e-15;
        _inner->dparm[9] = 1e-18;
        _inner->dparm[10] = 0.01;
        _inner->dparm[11] = 2e-8;
        _inner->dparm[21] = 2e-8;
        _inner->dparm[24] = 5e6;
        _inner->dparm[25] = 1.0;
        _inner->dparm[26] = 1.0;
        for (int i = 34; i < 63; i++)
        { _inner->dparm[i] = 0.0; }
        */

        if (fillStructure) {
            _inner->iparm[1] = 1;
            _inner->iparm[2] = 4;

            _inner->rowData.clear();
            _inner->rowData.reserve(_nrows + 1);

            // row data have to be always renumbered
            for (esint i = _data.K.nhalo; i <= _data.K.nrows; i++) {
                _inner->rowData.push_back(_data.K.rows[i] - _data.K.rows[_data.K.nhalo] + 1);
            }
            _inner->colIndices = _data.K.cols + _data.K.rows[_data.K.nhalo] - 1;
        }
        _inner->values = _data.K.vals + _data.K.rows[_data.K.nhalo] - 1;
    }
    _inner->rhsValues = _data.f[0].vals + _data.f[0].nhalo;

    eslog::solver("     - | PREPROCESSING                                                      %8.3f s | -\n", eslog::time() - start);
#endif
}

void WSMPSystemSolver::solve()
{
#ifdef HAVE_WSMP
    _inner->solution = _data.x[0].vals + _data.x[0].nhalo;
    double start = eslog::time();
    memcpy(_inner->solution, _inner->rhsValues, _nrows * sizeof(double));

    if (_inner->mtype == 1)
    {
        pwgsmp_(&_nrows, _inner->rowData.data(), _inner->colIndices,
            _inner->values, _inner->solution, &_inner->ldb,
            &_inner->nrhs, _inner->rmisc, _inner->iparm, _inner->dparm);
    } else {
        pwssmp_(&_nrows,_inner->rowData.data(), _inner->colIndices,
            _inner->values, _inner->diag, _inner->perm, _inner->invp,
            _inner->solution, &_inner->ldb, &_inner->nrhs, &_inner->aux,
            &_inner->naux, _inner->mrp, _inner->iparm, _inner->dparm);
    }

    _data.x.scatterToUpper();

    eslog::solver("     - | SOLVER TIME                                                        %8.3f s | -\n", eslog::time() - start);
    eslog::solver("     - --------------------------------------------------------------------------------- -\n");
#endif
}

WSMPSystemSolver::~WSMPSystemSolver()
{
#ifdef HAVE_WSMP
    delete[] _inner->rmisc;
    delete[] _inner->diag;
    delete[] _inner->mrp;
    delete[] _inner->invp;
    delete[] _inner->perm;

    delete _inner;
#endif
}

